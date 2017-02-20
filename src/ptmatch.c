#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mpfit.h"

#include "input.h"

// weighted distance between observed points and mapped reference points
int mapdist(int m, int n, double* p, double* d, double** dd, void* private)
{
    // problem description
    const int ni = n/5;
    const int nx = m/(ni - 1)/2;
    const double* x = private;
    
    // anchor point of image 0
    const double* a0 = p;
    
    // reference shear
    const double g1 = p[3];
    const double g2 = p[4];
    
    // number of deviates computed
    int k = 0;
    
    // try to map points from image 0 to image i
    for(int i = 1; i < ni; ++i)
    {
        // anchor point of image
        const double x0[2] = { p[5*i+0], p[5*i+1] };
        
        // matrix coefficients
        const double A = p[5*i+2];
        const double B = p[5*i+3];
        const double C = g2*A - g1*B;
        const double D = p[5*i+4];
        
        // transformation matrix
        const double T[4] = { 0.5*(D + A), 0.5*(B - C),
                              0.5*(B + C), 0.5*(D - A) };
        
        // map points, apply whitening transform and compute distance
        for(int j = 0; j < nx; ++j)
        {
            // reference point
            const double aj[2] = { x[5*nx*0+5*j+0], x[5*nx*0+5*j+1] };
            
            // observed point
            const double xj[2] = { x[5*nx*i+5*j+0], x[5*nx*i+5*j+1] };
            
            // whitening transform for observed point
            const double W[4] = { x[5*nx*i+5*j+2], x[5*nx*i+5*j+3],
                                                0, x[5*nx*i+5*j+4] };
            
            // offset between predicted and observed position
            const double delta[2] = {
                (xj[0] - x0[0]) - T[0]*(aj[0] - a0[0]) - T[1]*(aj[1] - a0[1]),
                (xj[1] - x0[1]) - T[2]*(aj[0] - a0[0]) - T[3]*(aj[1] - a0[1])
            };
            
            // compute uncorrelated deviates
            d[k+0] = W[0]*delta[0] + W[1]*delta[1];
            d[k+1] = W[2]*delta[0] + W[3]*delta[1];
            
            // compute derivatives if asked to
            if(dd)
            {
                // derivatives for reference shear
                if(dd[3])
                {
                    dd[3][k+0] = -0.5*B*(+ W[0]*(aj[1] - a0[1])
                                         - W[1]*(aj[0] - a0[0]));
                    dd[3][k+1] = -0.5*B*(+ W[2]*(aj[1] - a0[1])
                                         - W[3]*(aj[0] - a0[0]));
                }
                if(dd[4])
                {
                    dd[4][k+0] = +0.5*A*(+ W[0]*(aj[1] - a0[1])
                                         - W[1]*(aj[0] - a0[0]));
                    dd[4][k+1] = +0.5*A*(+ W[2]*(aj[1] - a0[1])
                                         - W[3]*(aj[0] - a0[0]));
                }
                
                // derivatives for anchor point
                if(dd[5*i+0])
                {
                    dd[5*i+0][k+0] = -W[0];
                    dd[5*i+0][k+1] = -W[2];
                }
                if(dd[5*i+1])
                {
                    dd[5*i+1][k+0] = -W[1];
                    dd[5*i+1][k+1] = -W[3];
                }
                
                // derivatives for a, b, d coefficients
                if(dd[5*i+2])
                {
                    dd[5*i+2][k+0] = -0.5*(+ W[0]*(aj[0] - a0[0])
                                           - W[0]*g2*(aj[1] - a0[1])
                                           - W[1]*(aj[1] - a0[1])
                                           + W[1]*g2*(aj[0] - a0[0]));
                    dd[5*i+2][k+1] = -0.5*(+ W[2]*(aj[0] - a0[0])
                                           - W[2]*g2*(aj[1] - a0[1])
                                           - W[3]*(aj[1] - a0[1])
                                           + W[3]*g2*(aj[0] - a0[0]));
                }
                if(dd[5*i+3])
                {
                    dd[5*i+3][k+0] = -0.5*(+ W[0]*(1 + g1)*(aj[1] - a0[1])
                                           + W[1]*(1 - g1)*(aj[0] - a0[0]));
                    dd[5*i+3][k+1] = -0.5*(+ W[2]*(1 + g1)*(aj[1] - a0[1])
                                           + W[3]*(1 - g1)*(aj[0] - a0[0]));
                }
                if(dd[5*i+4])
                {
                    dd[5*i+4][k+0] = -0.5*(+ W[0]*(aj[0] - a0[0])
                                           + W[1]*(aj[1] - a0[1]));
                    dd[5*i+4][k+1] = -0.5*(+ W[2]*(aj[0] - a0[0])
                                           + W[3]*(aj[1] - a0[1]));
                }
            }
            
            // done with deviates
            k += 2;
        }
    }
    
    // success if all deviates have been computed
    return k == m ? 0 : -1;
}

int main(int argc, char* argv[])
{
    // error indicator
    int err;
    
    // input
    char* f;
    char* o;
    char* m;
    char* a;
    int v, I, ND, DD;
    
    // number of images and points
    int ni, nx;
    double* x;
    
    // parameters
    int np;
    double* p;
    double* s;
    
    // optimiser data
    int       nd;
    mp_par*   par;
    mp_config cfg = {0};
    mp_result res = {0};
    
    
    /*********
     * input *
     *********/
    
    // default arguments
    f = o = m = a = NULL;
    v = I = ND = DD = 0;
    
    // parse arguments
    err = 0;
    for(int i = 1; i < argc && !err; ++i)
    {
        if(argv[i][0] == '-')
        {
            // flags
            for(char* c = &argv[i][1]; *c; ++c)
            {
                // increase verbosity
                if(*c == 'v')
                {
                    v += 1;
                }
                // decrease verbosity
                else if(*c == 'q')
                {
                    v -= 1;
                }
                // number of iterations
                else if(*c == 'I')
                {
                    if(!I)
                    {
                        if(i + 1 < argc)
                            I = atoi(argv[++i]);
                        else
                            err = 1;
                    }
                    else
                    {
                        err = 1;
                    }
                }
                // output file
                else if(*c == 'o')
                {
                    if(!o)
                    {
                        if(i + 1 < argc)
                            o = argv[++i];
                        else
                            err = 1;
                    }
                    else
                    {
                        err = 1;
                    }
                }
                // output matrix file
                else if(*c == 'm')
                {
                    if(!m)
                    {
                        if(i + 1 < argc)
                            m = argv[++i];
                        else
                            err = 1;
                    }
                    else
                    {
                        err = 1;
                    }
                }
                // output anchor file
                else if(*c == 'a')
                {
                    if(!a)
                    {
                        if(i + 1 < argc)
                            a = argv[++i];
                        else
                            err = 1;
                    }
                    else
                    {
                        err = 1;
                    }
                }
                // numerical derivatives (undocumented)
                else if(*c == 'N')
                {
                    ND = 1;
                }
                // debug derivatives (undocumented)
                else if(*c == 'D')
                {
                    DD = 1;
                }
                // unknown flag
                else
                {
                    err = 1;
                }
            }
        }
        else
        {
            // positional arguments
            if(!f)
                f = argv[i];
            else
                err = 1;
        }
    }
    
    // make sure input file was given
    if(!f)
        err = 1;
    
    // check for input errors
    if(err)
    {
        fprintf(stderr, "usage: ptmatch [-vq] [-I MAXITER] [-o OUTFILE] "
                "[-m MATFILE] [-a ANCFILE] PTSFILE\n");
        return EXIT_FAILURE;
    }
    
    // read points from input file
    read_points(f, &ni, &nx, &x);
    
    // make sure that enough images were given
    if(ni < 3)
    {
        fprintf(stderr, "%s: needs at least three images\n", f);
        return EXIT_FAILURE;
    }
    
    // make sure that enough points were given
    if(nx < 3)
    {
        fprintf(stderr, "%s: needs at least three points\n", f);
        return EXIT_FAILURE;
    }
    
    // compute points' whitening transform from input covariances
    for(int i = 0; i < ni; ++i)
    {
        for(int j = 0; j < nx; ++j)
        {
            // covariance matrix
            const double s1  = x[nx*5*i+5*j+2];
            const double s2  = x[nx*5*i+5*j+3];
            const double rho = x[nx*5*i+5*j+4];
            
            // check that given covariance is sane
            if(s1 < 0 || s2 < 0 || rho*rho >= 1)
            {
                fprintf(stderr, "%s: invalid covariance matrix for point %d "
                        "of image %d\n", f, j, i);
                return EXIT_FAILURE;
            }
            
            // Cholesky decomposition of inverse covariance matrix
            x[nx*5*i+5*j+2] = 1/(sqrt(1-rho*rho)*s1);
            x[nx*5*i+5*j+3] = -rho/(sqrt(1-rho*rho)*s2);
            x[nx*5*i+5*j+4] = 1/s2;
        }
    }
    
    
    /**************
     * initialise *
     **************/
    
    // total number of parameters, including fixed ones
    np = 5*ni;
    
    // create arrays for parameters and uncertainties
    p = malloc(np*sizeof(double));
    s = malloc(np*sizeof(double));
    if(!p || !s)
    {
        perror(NULL);
        return EXIT_FAILURE;
    }
    
    // compute a,b,c,d coefficients from given points
    for(int i = 1; i < ni; ++i)
    {
        // displacements from first three points in image 0
        const double u[4] = {
            x[5*1+0]-x[5*0+0], x[5*1+1]-x[5*0+1],
            x[5*2+0]-x[5*0+0], x[5*2+1]-x[5*0+1]
        };
        
        // displacements from first three points in image i
        const double v[4] = {
            x[nx*5*i+5*1+0]-x[nx*5*i+5*0+0], x[nx*5*i+5*1+1]-x[nx*5*i+5*0+1],
            x[nx*5*i+5*2+0]-x[nx*5*i+5*0+0], x[nx*5*i+5*2+1]-x[nx*5*i+5*0+1]
        };
        
        // transformation matrix from two sets of displacements
        const double T[4] = {
            (u[3]*v[0] - u[1]*v[2])/(u[0]*u[3] - u[1]*u[2]),
            (u[0]*v[2] - u[2]*v[0])/(u[0]*u[3] - u[1]*u[2]),
            (u[3]*v[1] - u[1]*v[3])/(u[0]*u[3] - u[1]*u[2]),
            (u[0]*v[3] - u[2]*v[1])/(u[0]*u[3] - u[1]*u[2])
        };
        
        // make sure transformation is valid
        if(!isfinite(T[0]) || !isfinite(T[1]) ||
           !isfinite(T[2]) || !isfinite(T[3]))
        {
            fprintf(stderr, "%s: points of image %d do not generate a "
                    "transformation matrix\n", f, i);
            return EXIT_FAILURE;
        }
        
        // store c,a,b,d coefficients for matrix
        p[5*i+1] = T[2] - T[1];
        p[5*i+2] = T[0] - T[3];
        p[5*i+3] = T[2] + T[1];
        p[5*i+4] = T[0] + T[3];
    }
    
    // convergence ratio for image 0 is fixed to unity
    p[2] = 1;
    
    // average initial shear from all pairs of multiple images
    p[3] = 0;
    p[4] = 0;
    
    // compute mean reference shear from a,b,c coefficients of images
    for(int i = 1, n = 0; i < ni; ++i)
    {
        // a, b, c coefficients of image i
        const double ci = p[5*i+1];
        const double ai = p[5*i+2];
        const double bi = p[5*i+3];
        
        // go through pairs of images i, j
        for(int j = i + 1; j < ni; ++j)
        {
            // a, b, c coefficients of image j
            const double cj = p[5*j+1];
            const double aj = p[5*j+2];
            const double bj = p[5*j+3];
            
            // compute reference shear from images i, j
            const double g1 = (ai*cj - aj*ci)/(bi*aj - bj*ai);
            const double g2 = (bi*cj - bj*ci)/(bi*aj - bj*ai);
            
            // update mean if shear from i, j is sane
            if(isfinite(g1) && isfinite(g2))
            {
                const double x1 = g1 - p[3];
                const double x2 = g2 - p[4];
                
                n += 1;
                p[3] += x1/n;
                p[4] += x2/n;
                s[3] += x1*(g1 - p[3]);
                s[4] += x2*(g2 - p[4]);
            }
        }
    }
    
    // output initial shear and a,b,d coefficients if very very verbose
    if(v >= 2)
    {
        printf("initial reference shear:\n");
        printf("  g1  % 10.4f  % 10.4f\n", p[3], sqrt(s[3]*2/ni/(ni-3)));
        printf("  g2  % 10.4f  % 10.4f\n", p[4], sqrt(s[4]*2/ni/(ni-3)));
        printf("initial a,b,d coefficients:\n");
        for(int i = 1; i < ni; ++i)
            printf("  %-#3.0f % 10.4f  % 10.4f  % 10.4f\n",
                   1.*i, p[5*i+2], p[5*i+3], p[5*i+4]);
    }
    
    // set anchor point to centroid of points for each image
    for(int i = 0; i < ni; ++i)
    {
        double c[2] = { 0, 0 };
        for(int j = 0; j < nx; ++j)
        {
            c[0] += x[nx*5*i+5*j+0];
            c[1] += x[nx*5*i+5*j+1];
        }
        p[5*i+0] = c[0]/nx;
        p[5*i+1] = c[1]/nx;
    }
    
    // output initial anchor points if very very verbose
    if(v >= 2)
    {
        printf("initial anchor points:\n");
        for(int i = 0; i < ni; ++i)
            printf("  %-#3.0f % 10.4f  % 10.4f\n", 1.*i, p[5*i+0], p[5*i+1]);
    }
    
    
    /************
     * minimise *
     ************/
    
    // number of data points
    nd = 2*nx*(ni - 1);
    
    // parameter configuration
    par = calloc(np, sizeof(mp_par));
    if(!par)
    {
        perror(NULL);
        return EXIT_FAILURE;
    }
    
    // fix reference image anchor point and convergence ratio
    par[0].fixed = 1;
    par[1].fixed = 1;
    par[2].fixed = 1;
    
    // set up other parameters
    for(int i = 3; i < np; ++i)
    {
        // numerical or analytical derivatives
        par[i].side = ND ? 0 : 3;
        
        // debug derivatives
        par[i].deriv_debug = DD;
    }
    
    // MPFIT configuration
    cfg.maxiter = I;
    cfg.nprint  = (v >= 2 ? 1 : 0);
    
    // set up results structure
    res.xerror = s;
    
    // minimise weighted distance between observed and mapped points
    err = mpfit(mapdist, nd, np, p, par, &cfg, x, &res);
    
    // check for errors
    if(err <= 0)
    {
        const char* msg = NULL;
        switch(err)
        {
            case MP_ERR_INPUT:
                msg = "General input parameter error";
                break;
            case MP_ERR_NAN:
                msg = "User function produced non-finite values";
                break;
            case MP_ERR_FUNC:
                msg = "No user function was supplied";
                break;
            case MP_ERR_NPOINTS:
                msg = "No user data points were supplied";
                break;
            case MP_ERR_NFREE:
                msg = "No free parameters";
                break;
            case MP_ERR_MEMORY:
                msg = "Memory allocation error";
                break;
            case MP_ERR_INITBOUNDS:
                msg = "Initial values inconsistent w constraint";
                break;
            case MP_ERR_BOUNDS:
                msg = "Initial constraints inconsistent";
                break;
            case MP_ERR_PARAM:
                msg = "General input parameter error";
                break;
            case MP_ERR_DOF:
                msg = "Not enough degrees of freedom";
                break;
        }
        if(msg)
            fprintf(stderr, "MPFIT error: %s \n", msg);
        else
            fprintf(stderr, "MPFIT error: %d \n", err);
        return EXIT_FAILURE;
    }
    
    // output convergence reason if verbose
    if(v >= 1)
    {
        const char* msg = NULL;
        switch(err)
        {
            case MP_OK_CHI:
                msg = "Convergence in chi-square value";
                break;
            case MP_OK_PAR:
                msg = "Convergence in parameter value";
                break;
            case MP_OK_BOTH:
                msg = "Both MP_OK_PAR and MP_OK_CHI hold";
                break;
            case MP_OK_DIR:
                msg = "Convergence in orthogonality";
                break;
            case MP_MAXITER:
                msg = "Maximum number of iterations reached";
                break;
            case MP_FTOL:
                msg = "ftol is too small; no further improvement";
                break;
            case MP_XTOL:
                msg = "xtol is too small; no further improvement";
                break;
            case MP_GTOL:
                msg = "gtol is too small; no further improvement";
                break;
        }
        if(msg)
            fprintf(stderr, "MPFIT success: %s \n", msg);
        else
            fprintf(stderr, "MPFIT success: %d \n", err);
    }
    
    // output results structure if very verbose
    if(v >= 2)
    {
        printf("MPFIT results:\n");
        printf("  bestnorm  = %16f\n", res.bestnorm);
        printf("  orignorm  = %16f\n", res.orignorm);
        printf("  niter     = %16d\n", res.niter);
        printf("  nfev      = %16d\n", res.nfev);
        printf("  status    = %16d\n", res.status);
        printf("  npar      = %16d\n", res.npar);
        printf("  nfree     = %16d\n", res.nfree);
        printf("  npegged   = %16d\n", res.npegged);
        printf("  nfunc     = %16d\n", res.nfunc);
    }
    
    
    /**********
     * output *
     **********/
    
    // convert g,a,b,d parameters to convergence ratio and shear
    for(int i = 1; i < ni; ++i)
    {
        // reference shear
        const double g1 = p[3];
        const double g2 = p[4];
        
        // a,b,c,d coefficients for image
        const double A = p[5*i+2];
        const double B = p[5*i+3];
        const double C = g2*A - g1*B;
        const double D = p[5*i+4];
        
        // numerator and denominator for f and g
        const double J = 0.5*(C*C + D*D - A*A - B*B);
        const double P = D*g1 - C*g2 + A;
        const double Q = C*g1 + D*g2 + B;
        const double R = A*g1 + B*g2 + D;
        
        // store initial f,g1,g2 for image
        p[5*i+2] = J/R;
        p[5*i+3] = P/R;
        p[5*i+4] = Q/R;
    }
    
    // print table of convergence ratios and shears
    if(v >= 0)
    {
        printf("        %10s  %10s\n", "ML", "sigma");
        for(int i = 0; i < ni; ++i)
        {
            printf("%-#3.0f  f  % 10.4f  % 10.4f\n", 1.*i, p[5*i+2], s[5*i+2]);
            printf("    g1  % 10.4f  % 10.4f\n", p[5*i+3], s[5*i+3]);
            printf("    g2  % 10.4f  % 10.4f\n", p[5*i+4], s[5*i+4]);
        }
    }
    
    // write convergence ratios and shears if asked to
    if(o)
    {
        // open matrix file for writing
        FILE* fp = fopen(o, "w");
        if(!fp)
        {
            perror(o);
            return EXIT_FAILURE;
        }
        
        // write convergence ratios and shears
        for(int i = 0; i < ni; ++i)
            fprintf(fp, "% 18.8f % 18.8f % 18.8f\n",
                    p[5*i+2], p[5*i+3], p[5*i+4]);
        
        // done
        fclose(fp);
    }
    
    // write transformation matrices if asked to
    if(m)
    {
        // open matrix file for writing
        FILE* fp = fopen(m, "w");
        if(!fp)
        {
            perror(m);
            return EXIT_FAILURE;
        }
        
        // write transformation matrices
        for(int i = 1; i < ni; ++i)
        {
            // convergence ratios and shears of images 0 and i
            const double g1 = p[3];
            const double g2 = p[4];
            const double F  = p[5*i+2];
            const double G1 = p[5*i+3];
            const double G2 = p[5*i+4];
        
            // matrix coefficients
            const double
            A = F*((1 - g1)*(1 + G1) - g2*G2)/(1 - G1*G1 - G2*G2),
            B = F*((1 + g1)*G2 - (1 + G1)*g2)/(1 - G1*G1 - G2*G2),
            C = F*((1 - g1)*G2 - (1 - G1)*g2)/(1 - G1*G1 - G2*G2),
            D = F*((1 + g1)*(1 - G1) - g2*G2)/(1 - G1*G1 - G2*G2);
            
            // write matrix coefficients
            fprintf(fp, "% 18.8f % 18.8f % 18.8f % 18.8f\n", A, B, C, D);
        }
        
        // done writing transforms
        fclose(fp);
    }
    
    // write anchor points if asked to
    if(a)
    {
        // open anchor file for writing
        FILE* fp = fopen(a, "w");
        if(!fp)
        {
            perror(m);
            return EXIT_FAILURE;
        }
        
        // write anchor points
        for(int i = 0; i < ni; ++i)
            fprintf(fp, "% 18.8f % 18.8f\n", p[5*i+0], p[5*i+1]);
        
        // done with anchor points file
        fclose(fp);
    }
    
    
    /************
     * cleaning *
     ************/
    
    free(x);
    free(p);
    free(s);
    free(par);
    
    
    /********
     * done *
     ********/
    
    return EXIT_SUCCESS;
}

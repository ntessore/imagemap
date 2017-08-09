#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "mpfit.h"
#include "mpis.h"
#include "input.h"

const char* USAGE =
"usage: ptmatch [-vqu] [-I MAXITER] [-o OUTFILE] [-m MATFILE]\n"
"               [-a ANCFILE] [-n NSAMPLE] [-s SAMFILE] PTSFILE";

// parameters during fit
enum {
    PAR_X,
    PAR_Y,
    PAR_A,
    PAR_B,
    PAR_C,
    PAR_D,
// ---
// last entry: number of parameters
// ---
    NP
};

// physical quantities
enum {
    PHYS_X,
    PHYS_Y,
    PHYS_J,
    PHYS_F,
    PHYS_G1,
    PHYS_G2
};

// weighted distance between observed points and mapped reference points
int mapdist(int m, int n, double* p, double* d, double** dd, void* private)
{
    // problem description
    const int ni = n/NP;
    const int nx = m/(ni - 1)/2;
    const double* x = private;
    
    // reference shear
    const double g01 = p[PHYS_G1];
    const double g02 = p[PHYS_G2];
    
    // number of deviates computed
    int k = 0;
    
    // try to map points from image 0 to image i
    for(int i = 1; i < ni; ++i)
    {
        // matrix coefficients
        const double ai = p[NP*i+PAR_A];
        const double bi = p[NP*i+PAR_B];
        const double ci = g02*ai - g01*bi;
        const double di = p[NP*i+PAR_D];
        
        // transformation matrix
        const double Ti[4] = { 0.5*(di + ai), 0.5*(bi - ci),
                               0.5*(bi + ci), 0.5*(di - ai) };
        
        // map points, apply whitening transform and compute distance
        for(int j = 0; j < nx; ++j)
        {
            // reference image point delta
            const double x0j[2] = { x[5*nx*0+5*j+0] - p[NP*0+PAR_X],
                                    x[5*nx*0+5*j+1] - p[NP*0+PAR_Y] };
            
            // multiple image point delta
            const double xij[2] = { x[5*nx*i+5*j+0] - p[NP*i+PAR_X],
                                    x[5*nx*i+5*j+1] - p[NP*i+PAR_Y] };
            
            // whitening transform for observed point
            const double Wij[4] = { x[5*nx*i+5*j+2], x[5*nx*i+5*j+3],
                                                  0, x[5*nx*i+5*j+4] };
            
            // offset between predicted and observed position
            const double Dij[2] = { xij[0] - Ti[0]*x0j[0] - Ti[1]*x0j[1],
                                    xij[1] - Ti[2]*x0j[0] - Ti[3]*x0j[1] };
            
            // compute uncorrelated deviates
            d[k+0] = Wij[0]*Dij[0] + Wij[1]*Dij[1];
            d[k+1] = Wij[2]*Dij[0] + Wij[3]*Dij[1];
            
            // compute derivatives if asked to
            if(dd)
            {
                // derivatives for reference shear
                if(dd[PHYS_G1])
                {
                    dd[PHYS_G1][k+0] = -0.5*bi*(Wij[0]*x0j[1] - Wij[1]*x0j[0]);
                    dd[PHYS_G1][k+1] = -0.5*bi*(Wij[2]*x0j[1] - Wij[3]*x0j[0]);
                }
                if(dd[PHYS_G2])
                {
                    dd[PHYS_G2][k+0] = +0.5*ai*(Wij[0]*x0j[1] - Wij[1]*x0j[0]);
                    dd[PHYS_G2][k+1] = +0.5*ai*(Wij[2]*x0j[1] - Wij[3]*x0j[0]);
                }
                
                // derivatives for anchor point
                if(dd[NP*i+PAR_X])
                {
                    dd[NP*i+0][k+0] = -Wij[0];
                    dd[NP*i+0][k+1] = -Wij[2];
                }
                if(dd[NP*i+PAR_Y])
                {
                    dd[NP*i+1][k+0] = -Wij[1];
                    dd[NP*i+1][k+1] = -Wij[3];
                }
                
                // derivatives for a, b, d coefficients
                if(dd[NP*i+PAR_A])
                {
                    dd[NP*i+PAR_A][k+0] = -0.5*(Wij[0]*x0j[0] - Wij[1]*x0j[1]
                                                - Wij[0]*g02*x0j[1]
                                                + Wij[1]*g02*x0j[0]);
                    dd[NP*i+PAR_A][k+1] = -0.5*(Wij[2]*x0j[0] - Wij[3]*x0j[1]
                                                - Wij[2]*g02*x0j[1]
                                                + Wij[3]*g02*x0j[0]);
                }
                if(dd[NP*i+PAR_B])
                {
                    dd[NP*i+PAR_B][k+0] = -0.5*(Wij[0]*(1 + g01)*x0j[1]
                                                + Wij[1]*(1 - g01)*x0j[0]);
                    dd[NP*i+PAR_B][k+1] = -0.5*(Wij[2]*(1 + g01)*x0j[1]
                                                + Wij[3]*(1 - g01)*x0j[0]);
                }
                if(dd[NP*i+PAR_D])
                {
                    dd[NP*i+PAR_D][k+0] = -0.5*(Wij[0]*x0j[0] + Wij[1]*x0j[1]);
                    dd[NP*i+PAR_D][k+1] = -0.5*(Wij[2]*x0j[0] + Wij[3]*x0j[1]);
                }
            }
            
            // done with deviates
            k += 2;
        }
    }
    
    // success if all deviates have been computed
    return k == m ? 0 : -1;
}

// convert parameters to physical quantities
void phys(int ni, double p[])
{
    // get reduced shear in reference image 0
    const double g01 = p[PHYS_G1];
    const double g02 = p[PHYS_G2];
    
    // magnification and convergence ratios in image 0 are unity
    p[NP*0+PHYS_J] = 1;
    p[NP*0+PHYS_F] = 1;
    
    // convert images to physical parameters
    for(int i = 1; i < ni; ++i)
    {
        // a,b,c,d coefficients for image i
        const double ai = p[NP*i+PAR_A];
        const double bi = p[NP*i+PAR_B];
        const double ci = g02*ai - g01*bi;
        const double di = p[NP*i+PAR_D];
        
        // determinant of relative magnification matrix
        const double J = 0.25*(ci*ci + di*di - ai*ai - bi*bi);
        
        // numerator and denominator for f and g
        const double P = di*g01 - ci*g02 + ai;
        const double Q = ci*g01 + di*g02 + bi;
        const double R = ai*g01 + bi*g02 + di;
        
        // compute j, f, g1, g2 for image i
        p[NP*i+PHYS_J]  = J;
        p[NP*i+PHYS_F]  = 2*J/R;
        p[NP*i+PHYS_G1] = P/R;
        p[NP*i+PHYS_G2] = Q/R;
    }
}

int main(int argc, char* argv[])
{
    // error indicator and message
    int err = 0;
    const char* msg = NULL;
    
    // input
    char* ptsfile;
    char* outfile;
    char* matfile;
    char* ancfile;
    char* samfile;
    int v, no_uncert, maxiter, nsample, ND, DD, seed;
    
    // number of images and points
    int ni, nx;
    double* x;
    
    // parameters
    int np;
    double* p;
    double* cov;
    
    // expectation mode
    int ns;
    double* s;
    double* m;
    double* e;
    
    // optimiser data
    int       nd;
    mp_par*   par;
    mp_config cfg = {0};
    mp_result res = {0};
    
    
    /*********
     * input *
     *********/
    
    // default arguments
    ptsfile = outfile = matfile = ancfile = samfile = NULL;
    v = no_uncert = maxiter = nsample = ND = DD = seed = 0;
    
    // parse arguments
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
                // ignore uncertainties
                else if(*c == 'u')
                {
                    if(!no_uncert)
                        no_uncert = 1;
                    else
                        err = 1;
                }
                // number of iterations
                else if(*c == 'I')
                {
                    if(!maxiter && i + 1 < argc)
                        maxiter = atoi(argv[++i]);
                    else
                        err = 1;
                }
                // output file
                else if(*c == 'o')
                {
                    if(!outfile && i + 1 < argc)
                        outfile = argv[++i];
                    else
                        err = 1;
                }
                // output matrix file
                else if(*c == 'm')
                {
                    if(!matfile && i + 1 < argc)
                        matfile = argv[++i];
                    else
                        err = 1;
                }
                // output anchor file
                else if(*c == 'a')
                {
                    if(!ancfile && i + 1 < argc)
                        ancfile = argv[++i];
                    else
                        err = 1;
                }
                // number of samples
                else if(*c == 'n')
                {
                    if(!nsample && i + 1 < argc)
                        nsample = atoi(argv[++i]);
                    else
                        err = 1;
                }
                // output samples file
                else if(*c == 's')
                {
                    if(!samfile && i + 1 < argc)
                        samfile = argv[++i];
                    else
                        err = 1;
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
                // seed (undocumented)
                else if(*c == 'S')
                {
                    if(!seed && i + 1 < argc)
                        seed = atoi(argv[++i]);
                    else
                        err = 1;
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
            if(!ptsfile)
                ptsfile = argv[i];
            else
                err = 1;
        }
    }
    
    // make sure input file was given
    if(!ptsfile)
        err = 1;
    
    // check for input errors
    if(err)
        goto err_usage;
    
    // read points from input file
    read_points(ptsfile, &ni, &nx, &x);
    
    // make sure that enough images were given
    if(ni < 3)
    {
        fprintf(stderr, "%s: needs at least three images\n", ptsfile);
        return EXIT_FAILURE;
    }
    
    // make sure that enough points were given
    if(nx < 3)
    {
        fprintf(stderr, "%s: needs at least three points\n", ptsfile);
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
                        "of image %d\n", ptsfile, j, i);
                return EXIT_FAILURE;
            }
            
            // check if uncertainties should be ignored
            if(no_uncert)
            {
                // default values
                x[nx*5*i+5*j+2] = 1;
                x[nx*5*i+5*j+3] = 0;
                x[nx*5*i+5*j+4] = 1;
            }
            else
            {
                // Cholesky decomposition of inverse covariance matrix
                x[nx*5*i+5*j+2] = 1/(sqrt(1-rho*rho)*s1);
                x[nx*5*i+5*j+3] = -rho/(sqrt(1-rho*rho)*s2);
                x[nx*5*i+5*j+4] = 1/s2;
            }
        }
    }
    
    
    /**************
     * initialise *
     **************/
    
    // total number of parameters, including fixed ones
    np = NP*ni;
    
    // create arrays for parameters and covariance matrix
    p = malloc(np*sizeof(double));
    cov = malloc(np*np*sizeof(double));
    if(!p || !cov)
        goto err_malloc;
    
    // set anchor point to centroid of points for all images
    for(int i = 0; i < ni; ++i)
    {
        p[NP*i+PAR_X] = 0;
        p[NP*i+PAR_Y] = 0;
        for(int j = 0; j < nx; ++j)
        {
            p[NP*i+PAR_X] += (x[nx*5*i+5*j+0] - p[NP*i+PAR_X])/(j + 1);
            p[NP*i+PAR_Y] += (x[nx*5*i+5*j+1] - p[NP*i+PAR_Y])/(j + 1);
        }
    }
    
    // compute a,b,c,d coefficients from given points
    for(int i = 1; i < ni; ++i)
    {
        // transformation matrix for image
        double T[4] = {0, 0, 0, 0};
        
        // number of transformations in mean
        int nt = 0;
        
        // go through points around centroid
        for(int j = 0, k = 1; j < nx; j = j + 1, k = (k + 1)%nx)
        {
            // offsets of reference points and anchor point in image 0
            const double u[4] = {
                x[nx*5*0+5*j+0] - p[NP*0+PAR_X],
                x[nx*5*0+5*j+1] - p[NP*0+PAR_Y],
                
                x[nx*5*0+5*k+0] - p[NP*0+PAR_X],
                x[nx*5*0+5*k+1] - p[NP*0+PAR_Y]
            };
            
            // offsets of reference points and anchor point in image i
            const double v[4] = {
                x[nx*5*i+5*j+0] - p[NP*i+PAR_X],
                x[nx*5*i+5*j+1] - p[NP*i+PAR_Y],
                
                x[nx*5*i+5*k+0] - p[NP*i+PAR_X],
                x[nx*5*i+5*k+1] - p[NP*i+PAR_Y]
            };
            
            // transformation matrix from two sets of displacements
            const double X[4] = {
                (u[3]*v[0] - u[1]*v[2])/(u[0]*u[3] - u[1]*u[2]),
                (u[0]*v[2] - u[2]*v[0])/(u[0]*u[3] - u[1]*u[2]),
                (u[3]*v[1] - u[1]*v[3])/(u[0]*u[3] - u[1]*u[2]),
                (u[0]*v[3] - u[2]*v[1])/(u[0]*u[3] - u[1]*u[2])
            };
            
            // add to mean if transformation is valid
            if(isfinite(X[0]) && isfinite(X[1]) &&
                isfinite(X[2]) && isfinite(X[3]))
            {
                nt += 1;
                T[0] += (X[0] - T[0])/nt;
                T[1] += (X[1] - T[1])/nt;
                T[2] += (X[2] - T[2])/nt;
                T[3] += (X[3] - T[3])/nt;
            }
        }
        
        // make sure a transformation matrix was found
        if(nt == 0)
        {
            fprintf(stderr, "%s: points of image %d do not generate a "
                    "transformation matrix\n", ptsfile, i);
            return EXIT_FAILURE;
        }
        
        // store coefficients a,b,c,d for matrix
        p[NP*i+PAR_A] = T[0] - T[3];
        p[NP*i+PAR_B] = T[2] + T[1];
        p[NP*i+PAR_C] = T[2] - T[1];
        p[NP*i+PAR_D] = T[0] + T[3];
    }
    
    // compute mean reference shear from a,b,c coefficients of images
    p[PHYS_G1] = 0;
    p[PHYS_G2] = 0;
    for(int i = 1, ng = 0; i < ni; ++i)
    {
        // a, b, c coefficients of image i
        const double ai = p[NP*i+PAR_A];
        const double bi = p[NP*i+PAR_B];
        const double ci = p[NP*i+PAR_C];
        
        // go through pairs of images i, j
        for(int j = i + 1; j < ni; ++j)
        {
            // a, b, c coefficients of image j
            const double aj = p[NP*j+PAR_A];
            const double bj = p[NP*j+PAR_B];
            const double cj = p[NP*j+PAR_C];
            
            // compute reference shear from images i, j
            const double g1 = (ai*cj - aj*ci)/(bi*aj - bj*ai);
            const double g2 = (bi*cj - bj*ci)/(bi*aj - bj*ai);
            
            // update mean if shear from i, j is sane
            if(isfinite(g1) && isfinite(g2))
            {
                ng += 1;
                p[PHYS_G1] += (g1 - p[PHYS_G1])/ng;
                p[PHYS_G2] += (g2 - p[PHYS_G2])/ng;
            }
        }
    }
    
    // output initial guess if very verbose
    if(v > 1)
    {
        double* q = malloc(np*sizeof(double));
        if(!q)
            goto err_malloc;
        for(int i = 0; i < np; ++i)
            q[i] = p[i];
        phys(ni, q);
        printf("initial values:\n");
        for(int i = 0; i < ni; ++i)
        {
            for(int j = 0; j < NP; ++j)
                printf("  % 10.4f", q[NP*i+j]);
            printf("\n");
        }
        free(q);
    }
    
    
    /************
     * minimise *
     ************/
    
    // number of data points
    nd = 2*nx*(ni - 1);
    
    // parameter configuration
    par = calloc(np, sizeof(mp_par));
    if(!par)
        goto err_malloc;
    
    // fix everything except reduced shear in reference image 0
    for(int j = 0; j < NP; ++j)
    {
        if(j == PHYS_G1 || j == PHYS_G2)
            par[j].fixed = 0;
        else
            par[j].fixed = 1;
    }
    
    // fix coefficient c in multiple images i
    for(int i = 1; i < ni; ++i)
        par[NP*i+PAR_C].fixed = 1;
    
    // set up parameters
    for(int i = 0; i < np; ++i)
    {
        // numerical or analytical derivatives
        par[i].side = ND ? 0 : 3;
        
        // debug derivatives
        par[i].deriv_debug = DD;
    }
    
    // MPFIT configuration
    cfg.maxiter = maxiter;
    cfg.nprint  = (v > 1 ? 1 : 0);
    
    // set up results structure
    res.covar = cov;
    
    // minimise weighted distance between observed and mapped points
    err = mpfit(mapdist, nd, np, p, par, &cfg, x, &res);
    
    // check for errors
    if(err <= 0)
        goto err_mpfit;
    
    // output convergence reason if verbose
    if(v > 0)
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
            printf("MPFIT success: %s \n", msg);
        else
            printf("MPFIT success: %d \n", err);
    }
    
    // output results structure if very verbose
    if(v > 1)
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
    
    
    /************
     * sampling *
     ************/
    
    // number of samples, use default if not specified
    ns = nsample ? nsample : 10000;
    
    // allocate space for samples
    s = malloc(ns*(np+2)*sizeof(double));
    if(!s)
        goto err_malloc;
    
    // seed random number generator
    srand(seed ? seed : time(0));
    
    // status output if verbose
    if(v > 0)
        printf("expectation mode: %d samples\n", ns);
    
    // draw samples
    err = mpis(mapdist, nd, np, p, cov, x, ns, s);
    if(err <= 0)
        goto err_mpfit;
    
    // output effective number of samples if verbose
    if(v > 0)
        printf("effective number of samples: %d\n",
               (int)mpis_neff(np, ns, s));
    
    
    /***********
     * results *
     ***********/
    
    // convert maximum-likelihood parameters to f and g
    phys(ni, p);
    
    // allocate space for sampling results
    m = malloc(np*sizeof(double));
    e = malloc(np*sizeof(double));
    if(!m || !e)
        goto err_malloc;
    
    // convert sampled parameters to f and g
    for(int i = 0; i < ns; ++i)
        phys(ni, &s[i*(np+2)]);
    
    // compute mean and error of lens quantities from samples
    mpis_stat(np, ns, s, m, e);
    
    
    /**********
     * output *
     **********/
    
    // print table of convergence ratios and shear if not quiet
    if(v > -1)
    {
        const char* labels[] = { "x", "y", "j", "f", "g1", "g2" };
        
        printf("+--------+------------+------------+------------+\n");
        printf("| %-6s | %10s | %10s | %10s |\n",
               "value", "ML", "mean", "sigma");
        printf("+--------+------------+------------+------------+\n");
        
        for(int i = 0; i < ni; ++i)
        {
            for(int j = 2; j < NP; ++j)
                printf("| %2s_%-3d | % 10.4f | % 10.4f | % 10.4f |\n",
                       labels[j], i, p[NP*i+j], m[NP*i+j], e[NP*i+j]);
            printf("+--------+------------+------------+------------+\n");
        }
    }
    
    // write convergence ratios and shears if asked to
    if(outfile)
    {
        // open matrix file for writing
        FILE* fp = fopen(outfile, "w");
        if(!fp)
        {
            msg = outfile;
            goto err_file;
        }
        
        // write convergence ratios and shears
        for(int i = 0; i < ni; ++i)
            fprintf(fp, "% 18.8f % 18.8f % 18.8f\n",
                    p[NP*i+PHYS_F], p[NP*i+PHYS_G1], p[NP*i+PHYS_G2]);
        
        // done
        fclose(fp);
    }
    
    // write transformation matrices if asked to
    if(matfile)
    {
        // open matrix file for writing
        FILE* fp = fopen(matfile, "w");
        if(!fp)
        {
            msg = matfile;
            goto err_file;
        }
        
        // write transformation matrices
        for(int i = 1; i < ni; ++i)
        {
            // convergence ratios and shears of images 0 and i
            const double g1 = p[3];
            const double g2 = p[4];
            const double F  = p[NP*i+PHYS_F];
            const double G1 = p[NP*i+PHYS_G1];
            const double G2 = p[NP*i+PHYS_G2];
            
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
    if(ancfile)
    {
        // open anchor file for writing
        FILE* fp = fopen(ancfile, "w");
        if(!fp)
        {
            msg = ancfile;
            goto err_file;
        }
        
        // write anchor points
        for(int i = 0; i < ni; ++i)
            fprintf(fp, "% 18.8f % 18.8f\n", p[NP*i+PHYS_X], p[NP*i+PHYS_Y]);
        
        // done with anchor points file
        fclose(fp);
    }
    
    // write samples if asked to
    if(samfile)
    {
        FILE* fp = fopen(samfile, "w");
        if(!fp)
        {
            msg = samfile;
            goto err_file;
        }
        for(int i = 0; i < ns; ++i)
        {
            fprintf(fp, "  %28.18E", s[i*(np+2)+np+1]);
            fprintf(fp, "  %28.18E", -2*s[i*(np+2)+np+0]);
            for(int j = 0; j < np; ++j)
                fprintf(fp, "  %28.18E", s[i*(np+2)+j]);
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    
    
    /************
     * cleaning *
     ************/
    
    free(x);
    free(p);
    free(cov);
    free(par);
    free(s);
    free(m);
    free(e);
    
    
    /********
     * done *
     ********/
    
    return EXIT_SUCCESS;
    
    
    /**********
     * errors *
     **********/
    
    // print usage message, might not be an error
err_usage:
    fprintf(err ? stderr : stdout, "%s\n", USAGE);
    if(err && msg)
        fprintf(stderr, "\nerror: %s\n", msg);
    return err ? EXIT_FAILURE : EXIT_SUCCESS;
    
    // memory allocation error
err_malloc:
    perror(NULL);
    return EXIT_FAILURE;
    
    // file error
err_file:
    perror(msg);
    return EXIT_FAILURE;
    
    // handle CMPFIT error codes
err_mpfit:
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
    
    // print error message or code
    if(msg)
        fprintf(stderr, "MPFIT error: %s\n", msg);
    else
        fprintf(stderr, "MPFIT error: %d\n", err);
    
    return EXIT_FAILURE;
}

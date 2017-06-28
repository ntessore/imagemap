#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "mpfit.h"
#include "input.h"

const char* USAGE =
"usage: ptmatch [-vqux] [-I MAXITER] [-o OUTFILE] [-m MATFILE]\n"
"               [-a ANCFILE] [-n NSAMPLE] [-s SAMFILE] PTSFILE";

// epsilon for stabilising the Cholesky decomposition
#ifndef CHOL_EPS
#define CHOL_EPS 1E-100
#endif

// Cholesky decomposition
void chol_eps(int n, const double A[], double eps, double L[])
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < i + 1; ++j)
        {
            double s = 0;
            for(int k = 0; k < j; ++k)
                s += L[i*n+k]*L[j*n+k];
            if(j < i)
                L[i*n+j] = (A[i*n+j] - s)/L[j*n+j];
            else
                L[i*n+j] = sqrt(A[i*n+i] + eps - s);
        }
        for(int j = i + 1; j < n; ++j)
            L[i*n+j] = 0;
    }
}

// generate array of random normal variates using the polar method
void randnv(int c, double v[])
{
    double x[2];
    double s = 0;
    
    for(int i = 0; i < c; ++i)
    {
        if(i%2 == 0)
        {
            do
            {
                x[0] = 2.*rand()/(1. + RAND_MAX) - 1;
                x[1] = 2.*rand()/(1. + RAND_MAX) - 1;
                s = x[0]*x[0] + x[1]*x[1];
            }
            while(s >= 1);
        }
        v[i] = x[i%2]*sqrt(-2*log(s)/s);
    }
}

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

// convert g,a,b,d parameters to convergence ratio f and reduced shear g
void ptofg(int ni, double p[])
{
    // reference shear
    const double g1 = p[3];
    const double g2 = p[4];
    
    for(int i = 1; i < ni; ++i)
    {
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
        
        // compute f, g1, g2 for image
        p[5*i+2] = J/R;
        p[5*i+3] = P/R;
        p[5*i+4] = Q/R;
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
    int v, xmode, no_uncert, maxiter, nsample, ND, DD, seed;
    
    // number of images and points
    int ni, nx;
    double* x;
    
    // parameters
    int np;
    double* p;
    double* cov;
    
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
    v = xmode = no_uncert = maxiter = nsample = ND = DD = seed = 0;
    
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
                // expectation mode
                else if(*c == 'x')
                {
                    if(!xmode)
                        xmode = 1;
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
    
    // sample output requires expectation mode
    if(samfile && !xmode)
    {
        msg = "sample output (-s) requires expectation mode (-x)";
        err = 1;
    }
    
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
    np = 5*ni;
    
    // create arrays for parameters and covariance matrix
    p = malloc(np*sizeof(double));
    cov = malloc(np*np*sizeof(double));
    if(!p || !cov)
        goto err_malloc;
    
    // set anchor point to centroid of points for image 0
    p[0] = p[1] = 0;
    for(int j = 0; j < nx; ++j)
    {
        p[0] += (x[5*j+0] - p[0])/(j + 1);
        p[1] += (x[5*j+1] - p[1])/(j + 1);
    }
    
    // compute a,b,c,d coefficients from given points
    for(int i = 1; i < ni; ++i)
    {
        // transformation matrix for image
        double T[4] = {0, 0, 0, 0};
        
        // number of transformations in mean
        int nt = 0;
        
        // centroid of points for image i
        double q[2] = {0, 0};
        for(int j = 0; j < nx; ++j)
        {
            q[0] += (x[nx*5*i+5*j+0] - q[0])/(j + 1);
            q[1] += (x[nx*5*i+5*j+1] - q[1])/(j + 1);
        }
        
        // go through points around centroid
        for(int j = 0, k = 1; j < nx; j = j + 1, k = (k + 1)%nx)
        {
            // offsets of reference points and anchor point in image 0
            const double u[4] = {
                x[nx*5*0+5*j+0] - p[0], x[nx*5*0+5*j+1] - p[1],
                x[nx*5*0+5*k+0] - p[0], x[nx*5*0+5*k+1] - p[1]
            };
            
            // offsets of reference points and anchor point in image i
            const double v[4] = {
                x[nx*5*i+5*j+0] - q[0], x[nx*5*i+5*j+1] - q[1],
                x[nx*5*i+5*k+0] - q[0], x[nx*5*i+5*k+1] - q[1]
            };
            
            // transformation matrix from two sets of displacements
            const double
            A = (u[3]*v[0] - u[1]*v[2])/(u[0]*u[3] - u[1]*u[2]),
            B = (u[0]*v[2] - u[2]*v[0])/(u[0]*u[3] - u[1]*u[2]),
            C = (u[3]*v[1] - u[1]*v[3])/(u[0]*u[3] - u[1]*u[2]),
            D = (u[0]*v[3] - u[2]*v[1])/(u[0]*u[3] - u[1]*u[2]);
            
            // add to mean if transformation is valid
            if(isfinite(A) && isfinite(B) && isfinite(C) && isfinite(D))
            {
                nt += 1;
                T[0] += (A - T[0])/nt;
                T[1] += (B - T[1])/nt;
                T[2] += (C - T[2])/nt;
                T[3] += (D - T[3])/nt;
            }
        }
        
        // make sure a transformation matrix was found
        if(nt == 0)
        {
            fprintf(stderr, "%s: points of image %d do not generate a "
                    "transformation matrix\n", ptsfile, i);
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
    
    // compute mean reference shear from a,b,c coefficients of images
    p[3] = p[4] = 0;
    for(int i = 1, ng = 0; i < ni; ++i)
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
                ng += 1;
                p[3] += (g1 - p[3])/ng;
                p[4] += (g2 - p[4])/ng;
            }
        }
    }
    
    // set anchor point to centroid of points for each multiple image
    for(int i = 1; i < ni; ++i)
    {
        p[5*i+0] = p[5*i+1] = 0;
        for(int j = 0; j < nx; ++j)
        {
            p[5*i+0] += (x[nx*5*i+5*j+0] - p[5*i+0])/(j + 1);
            p[5*i+1] += (x[nx*5*i+5*j+1] - p[5*i+1])/(j + 1);
        }
    }
    
    // output initial parameter guess if very verbose
    if(v > 1)
    {
        printf("initial parameters:\n");
        for(int i = 0; i < ni; ++i)
        {
            for(int j = 0; j < 5; ++j)
                printf("  % 10.4f", p[5*i+j]);
            printf("\n");
        }
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
    
    
    /********************
     * expectation mode *
     ********************/
    
    // compute expectation values when `-x` is given
    if(xmode)
    {
        // number of samples
        int ns;
        
        // sample, cumulative weight, cumulative square weight
        double w, W, W2;
        
        // sample log-likelihood
        double loglike;
        
        // array for mean parameter values
        double* m;
        
        // arrays for sampling (random normal variate, parameters, deviates)
        double* V;
        double* P;
        double* D;
        
        // sample file
        FILE* fp;
        
        // seed random number generator
        srand(seed ? seed : time(0));
        
        // number of samples, default is 10000
        ns = nsample ? nsample : 10000;
        
        // allocate zeroed array for mean parameter values
        m = calloc(np, sizeof(double));
        if(!m)
            goto err_malloc;
        
        // allocate arrays for sampling
        V = malloc(np*sizeof(double));
        P = malloc(np*sizeof(double));
        D = malloc(nd*sizeof(double));
        if(!V || !P || !D)
            goto err_malloc;
        
        // open sample file for writing if asked to
        if(samfile)
        {
            fp = fopen(samfile, "w");
            if(!fp)
            {
                msg = samfile;
                goto err_file;
            }
        }
        else
        {
            // no output
            fp = NULL;
        }
        
        // status output if verbose
        if(v > 0)
            printf("expectation mode: %d samples\n", ns);
        
        // compute square root of covariance matrix (in place)
        chol_eps(np, cov, CHOL_EPS, cov);
        
        // zero cumulative weight
        W = W2 = 0;
        
        // draw samples with given mean and variance
        for(int i = 0; i < ns; ++i)
        {
            // draw uncorrelated random normal variates
            randnv(np, V);
            
            // compute log-probability of draw, ignore fixed parameters
            w = 0;
            for(int j = 0; j < np; ++j)
                if(!par[j].fixed)
                    w += -0.5*V[j]*V[j];
            
            // transform using mean and sqrt(covar)
            for(int j = 0; j < np; ++j)
            {
                P[j] = p[j];
                for(int k = 0; k < np; ++k)
                    P[j] += cov[j*np+k]*V[k];
            }
            
            // compute deviates
            err = mapdist(nd, np, P, D, NULL, x);
            
            // check for user function error
            if(err < 0)
                goto err_mpfit;
            
            // compute log-likelihood of parameters
            loglike = 0;
            for(int j = 0; j < nd; ++j)
                loglike += -0.5*D[j]*D[j];
            
            // compute sample weight
            w = exp(loglike - w);
            
            // add to cumulative weight
            W += w;
            W2 += w*w;
            
            // convert sample parameters to f and g
            ptofg(ni, P);
            
            // add to weighted mean
            for(int j = 0; j < np; ++j)
                m[j] += (w/W)*(P[j] - m[j]);
            
            // write sample if asked to
            if(fp)
            {
                fprintf(fp, "%28.18E", w);
                fprintf(fp, "%28.18E", loglike);
                for(int j = 0; j < np; ++j)
                    fprintf(fp, "%28.18E", P[j]);
                fprintf(fp, "\n");
            }
        }
        
        // output effective number of samples if verbose
        if(v > 0)
            printf("effective number of samples: %d\n", (int)(W*W/W2));
        
        // close sample file if open
        if(fp)
            fclose(fp);
        
        // weighted means replace parameter values
        free(p);
        p = m;
        
        // done with sampling
        free(V);
        free(P);
        free(D);
    }
    else
    {
        // convert best-fit parameters to f and g
        ptofg(ni, p);
    }
    
    
    /**********
     * output *
     **********/
    
    // print table of convergence ratios and shear if not quiet
    if(v > -1)
    {
        for(int i = 0; i < ni; ++i)
            printf("% 18.8f % 18.8f % 18.8f\n", p[5*i+2], p[5*i+3], p[5*i+4]);
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
                    p[5*i+2], p[5*i+3], p[5*i+4]);
        
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
            fprintf(fp, "% 18.8f % 18.8f\n", p[5*i+0], p[5*i+1]);
        
        // done with anchor points file
        fclose(fp);
    }
    
    
    /************
     * cleaning *
     ************/
    
    free(x);
    free(p);
    free(cov);
    free(par);
    
    
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

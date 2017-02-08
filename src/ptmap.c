#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "newuoa.h"

#include "input.h"

// initial and final radius of trust region
static const double RBEG = 1e0;
static const double REND = 1e-10;

// maximum number of function evaluation
static const int FLIM = 1000000;

// context for optimiser
typedef struct
{
    int ni;
    int nx;
    double* x;
    double* a0;
} context;

// weighted distance between observed points and mapped reference points
double mapdist(const long n, const double* p, void* ctx_)
{
    // get context structure
    const context* ctx = ctx_;
    
    // number of images and points
    const int ni = ctx->ni;
    const int nx = ctx->nx;
    const double* x = ctx->x;
    
    // anchor point of image 0
    const double* a0 = ctx->a0;
    
    // reference shear
    const double g1 = p[0];
    const double g2 = p[1];
    
    // total distance between observed and mapped points
    double dist = 0;
    
    // try to map points from image 0 to image i
    for(int i = 1; i < ni; ++i)
    {
        // anchor point of image
        const double x0[2] = { p[5*i-3], p[5*i-2] };
        
        // matrix coefficients
        const double a = p[5*i-1];
        const double b = p[5*i+0];
        const double c = g2*a - g1*b;
        const double d = p[5*i+1];
        
        // transformation matrix
        const double T[4] = { 0.5*(d + a), 0.5*(b - c),
                              0.5*(b + c), 0.5*(d - a) };
        
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
            const double D[2] = {
                (xj[0] - x0[0]) - T[0]*(aj[0] - a0[0]) - T[1]*(aj[1] - a0[1]),
                (xj[1] - x0[1]) - T[2]*(aj[0] - a0[0]) - T[3]*(aj[1] - a0[1])
            };
            
            // compute uncorrelated deviates
            const double d[2] = {
                W[0]*D[0] + W[1]*D[1],
                W[2]*D[0] + W[3]*D[1]
            };
            
            // add to total distance
            dist += d[0]*d[0] + d[1]*d[1];
        }
    }
    
    // return the total distance
    return dist;
}

int main(int argc, char* argv[])
{
    // error indicator
    int err;
    
    // input
    char* f;
    char* o;
    char* m;
    int v;
    
    // number of images and points
    int ni;
    int nx;
    double* x;
    
    // parameters
    int np;
    double* p;
    
    // optimiser data
    int n, npt;
    context ctx;
    double* tmp;
    
    
    /*********
     * input *
     *********/
    
    // default arguments
    f = o = m = NULL;
    v = 0;
    
    // parse arguments
    err = 0;
    for(int i = 1; i < argc && !err; ++i)
    {
        if(argv[i][0] == '-')
        {
            // flags
            for(char* a = &argv[i][1]; *a; ++a)
            {
                // increase verbosity
                if(*a == 'v')
                {
                    v += 1;
                }
                // decrease verbosity
                else if(*a == 'q')
                {
                    v -= 1;
                }
                // output file
                else if(*a == 'o')
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
                else if(*a == 'm')
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
        fprintf(stderr, "usage: ptmap [-vq] [-o OUTFILE] [-m MATFILE] FILE\n");
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
    
    // create array for parameters
    p = malloc(np*sizeof(double));
    if(!p)
    {
        perror(NULL);
        return EXIT_FAILURE;
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
    
    // compute a,b,d coefficients from given points
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
        
        // store a,b,d coefficients for matrix
        p[5*i+2] = T[0] - T[3];
        p[5*i+3] = T[2] + T[1];
        p[5*i+4] = T[0] + T[3];
        
        // store c coefficient for images 1 and 2 only
        if(i < 3)
            p[2+i] = T[2] - T[1];
    }
    
    // convergence ratio for image 0 is unity by definition
    p[2] = 1;
    
    // compute shear for image 0 from a,b,c coefficients of images 1 and 2
    {
        // stored a, b, c coefficients of images 1 and 2
        const double a1 = p[5*1+2];
        const double b1 = p[5*1+3];
        const double c1 = p[2+1];
        const double a2 = p[5*2+2];
        const double b2 = p[5*2+3];
        const double c2 = p[2+2];
        
        // compute g for reference image 0
        p[3] = (a1*c2 - a2*c1)/(b1*a2 - b2*a1);
        p[4] = (b1*c2 - b2*c1)/(b1*a2 - b2*a1);
        
        // make sure that initial shear is sane
        if(!isfinite(p[3]) || !isfinite(p[4]))
        {
            fprintf(stderr, "%s: images 1 and 2 do not generate a valid "
                    "shear\n", f);
            return EXIT_FAILURE;
        }
    }
    
    
    /************
     * minimise *
     ************/
    
    // dimensionality of problem and interpolation
    n = np - 3;
    npt = 2*n + 1;
    
    // set up data structure
    ctx.ni = ni;
    ctx.nx = nx;
    ctx.x  = x;
    ctx.a0 = p;
    
    // allocate scratch space
    tmp = malloc(((npt+13)*(npt+n)+3*n*(n+3)/2)*sizeof(double));
    if(!tmp)
    {
        perror(NULL);
        return EXIT_FAILURE;
    }
    
    // minimise weighted distance between observed and mapped points
    err = newuoa(n, npt, mapdist, &ctx, &p[3], RBEG, REND, v, FLIM, tmp);
    
    // free scratch space
    free(tmp);
    
    // make sure call was successful
    if(err != NEWUOA_SUCCESS)
    {
        fprintf(stderr, "error in NEWUOA: %s\n", newuoa_reason(err));
        return EXIT_FAILURE;
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
        const double a = p[5*i+2];
        const double b = p[5*i+3];
        const double c = g2*a - g1*b;
        const double d = p[5*i+4];
        
        // numerator and denominator for f and g
        const double J = 0.5*(c*c + d*d - a*a - b*b);
        const double P = d*g1 - c*g2 + a;
        const double Q = c*g1 + d*g2 + b;
        const double R = a*g1 + b*g2 + d;
        
        // store initial f,g1,g2 for image
        p[5*i+2] = J/R;
        p[5*i+3] = P/R;
        p[5*i+4] = Q/R;
    }
    
    // print table of convergence ratios and shears
    if(v >= 0)
    {
        printf("%4s  %8s  %8s  %8s\n", "i", "f", "g_1", "g_2");
        for(int i = 0; i < ni; ++i)
            printf("%4d  %8.4f  %8.4f  %8.4f\n",
                   i, p[5*i+2], p[5*i+3], p[5*i+4]);
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
    
    
    /************
     * cleaning *
     ************/
    
    free(x);
    free(p);
    
    
    /********
     * done *
     ********/
    
    return EXIT_SUCCESS;
}

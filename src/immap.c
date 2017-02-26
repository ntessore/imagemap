#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fitsio.h"

#include "input.h"

const char* USAGE =
"usage: immap [-0] [-p PTSFILE] [-a ANCFILE] [-d PAD] MATFILE IMFITS OUTFITS";

// subsampling
#ifndef NSUB
#define NSUB 32
#endif

// value for null pixels
static const double NULVAL = DOUBLENULLVALUE;

int main(int argc, char* argv[])
{
    // error flag and message
    int err = 0;
    const char* msg = NULL;
    
    // input
    char* matfile;
    char* imfile;
    char* outfile;
    char* ptsfile;
    char* ancfile;
    int use_null;
    double pad;
    
    // matrix table
    int ni, nr, nc;
    double* mat;
    
    // points
    int nx;
    double* x;
    double* a;
    
    // FITS file
    fitsfile* fp;
    int na;
    long n[2], npix, fpix[2] = { 1, 1 };
    void* nulv;
    char extname[10];
    double* pix;
    
    // mapped image bounding box
    double A[2], B[2];
    
    // default arguments
    matfile = imfile = outfile = ptsfile = ancfile = NULL;
    use_null = 0;
    pad = 0.;
    
    // parse arguments
    err = 0;
    for(int i = 1; i < argc && !err; ++i)
    {
        if(argv[i][0] == '-')
        {
            // flags
            for(char* c = &argv[i][1]; *c; ++c)
            {
                // write null values instead of zeros
                if(*c == '0')
                {
                    if(!use_null)
                        use_null = 1;
                    else
                        err = 1;
                }
                // points file
                else if(*c == 'p')
                {
                    if(!ptsfile && i + 1 < argc)
                        ptsfile = argv[++i];
                    else
                        err = 1;
                }
                // anchor points file
                else if(*c == 'a')
                {
                    if(!ancfile && i + 1 < argc)
                        ancfile = argv[++i];
                    else
                        err = 1;
                }
                // mapped image padding
                else if(*c == 'd')
                {
                    if(!pad && i + 1 < argc)
                        pad = atof(argv[++i]);
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
            if(!matfile)
                matfile = argv[i];
            else if(!imfile)
                imfile = argv[i];
            else if(!outfile)
                outfile = argv[i];
            else
                err = 1;
        }
    }
    
    // make sure input files were given
    if(!matfile || !imfile || !outfile)
        err = 1;
    
    // anchor points need points
    if(ancfile && !ptsfile)
    {
        msg = "anchor points (-a) can only be given with points (-p)";
        err = 1;
    }
    
    // set default padding to 50%
    if(!pad)
        pad = -0.5;
    
    // check for input errors
    if(err)
        goto err_usage;
    
    // set null value if not writing zeros
    if(use_null)
        nulv = (void*)&NULVAL;
    else
        nulv = NULL;
    
    // read input file
    read_table(matfile, &nr, &nc, &mat);
    
    // make sure format is correct
    if(nc != 4)
    {
        msg = "file needs four columns (matrix entries)";
        goto err_matfile;
    }
    
    // number of images
    ni = nr + 1;
    
    // read points if given
    if(ptsfile)
    {
        // read points from file
        read_points(ptsfile, &nr, &nx, &x);
        
        // make sure number of images match
        if(nr != ni)
        {
            msg = "number of images in point file and matrix file must match";
            goto err_ptsfile;
        }
    }
    else
    {
        // no points given
        nx = 0;
        x = NULL;
    }
    
    // read anchor points if given
    if(ancfile)
    {
        // read file as table
        read_table(ancfile, &nr, &nc, &a);
        
        // check format
        if(nc != 2)
        {
            msg = "anchor points file needs two columns (x, y)";
            goto err_ancfile;
        }
        
        // check number of images
        if(nr != ni)
        {
            msg = "number of anchor points must match number of images";
            goto err_ancfile;
        }
    }
    else
    {
        // array for anchor points
        a = malloc(2*ni*sizeof(double));
        if(!a)
            goto err_malloc;
        
        // set anchor points
        for(int i = 0; i < ni; ++i)
        {
            // set to zero
            a[2*i+0] = 0;
            a[2*i+1] = 0;
            
            // if points are given, compute centroid
            if(nx > 0)
            {
                for(int j = 0; j < nx; ++j)
                {
                    a[2*i+0] += x[5*nx*i+5*j+0];
                    a[2*i+1] += x[5*nx*i+5*j+1];
                }
                a[2*i+0] /= nx;
                a[2*i+1] /= nx;
            }
        }
    }
    
    // read input FITS
    fits_open_image(&fp, imfile, READONLY, &err);
    fits_get_img_dim(fp, &na, &err);
    fits_get_img_size(fp, 2, n, &err);
    
    // make sure the image is an image
    if(!err && na != 2)
        err = BAD_DIMEN;
    
    // check for errors
    if(err)
        goto err_imfits;
    
    // allocate pixel array
    npix = n[0]*n[1];
    pix = malloc(npix*sizeof(double));
    if(!pix)
        goto err_malloc;
    
    // read pixels into array and close input file
    fits_read_pix(fp, TDOUBLE, fpix, npix, NULL, pix, NULL, &err);
    fits_close_file(fp, &err);
    
    // check for errors
    if(err)
        goto err_imfits;
    
    // create output FITS
    fits_create_file(&fp, outfile, &err);
    
    // check for errors
    if(err)
        goto err_outfits;
    
    // transform and write images
    for(int i = 0; i < ni; ++i)
    {
        // matrix for image
        double T[4];
        
        // no transformation for first image
        if(i == 0)
        {
            T[0] = 1;
            T[1] = 0;
            T[2] = 0;
            T[3] = 1;
        }
        else
        {
            T[0] = mat[4*i - 4];
            T[1] = mat[4*i - 3];
            T[2] = mat[4*i - 2];
            T[3] = mat[4*i - 1];
        };
        
        // determinant of matrix
        const double det = fabs(T[0]*T[3] - T[1]*T[2]);
        
        // get mapped image bounding box
        if(nx > 0)
        {
            // get bounding box from provided points
            for(int j = 0; j < nx; ++j)
            {
                // get point
                const double p[2] = {
                    x[5*nx*i + 5*j + 0],
                    x[5*nx*i + 5*j + 1]
                };
                
                // expand bounding box if necessary
                if(j == 0 || A[0] > p[0])
                    A[0] = p[0];
                if(j == 0 || A[1] > p[1])
                    A[1] = p[1];
                if(j == 0 || B[0] < p[0])
                    B[0] = p[0];
                if(j == 0 || B[1] < p[1])
                    B[1] = p[1];
            }
            
            // pad image
            if(pad > 0)
            {
                // positive number: pixel padding
                A[0] -= pad;
                A[1] -= pad;
                B[0] += pad;
                B[1] += pad;
            }
            else
            {
                // negative number: percent padding
                const double P[2] = {
                    pad*(A[0] - B[0]),
                    pad*(A[1] - B[1])
                };
                A[0] -= P[0];
                A[1] -= P[1];
                B[0] += P[0];
                B[1] += P[1];
            }
        }
        else
        {
            // map image corners
            const double C[4][2] = {
                { T[0]*1    + T[1]*1   , T[2]*1    + T[3]*1    },
                { T[0]*n[0] + T[1]*1   , T[2]*n[0] + T[3]*1    },
                { T[0]*1    + T[1]*n[1], T[2]*1    + T[3]*n[1] },
                { T[0]*n[0] + T[1]*n[1], T[2]*n[0] + T[3]*n[1] }
            };
        
            // get mapped image bounding box
            A[0] = fmin(C[0][0], fmin(C[1][0], fmin(C[2][0], C[3][0])));
            A[1] = fmin(C[0][1], fmin(C[1][1], fmin(C[2][1], C[3][1])));
            B[0] = fmax(C[0][0], fmax(C[1][0], fmax(C[2][0], C[3][0])));
            B[1] = fmax(C[0][1], fmax(C[1][1], fmax(C[2][1], C[3][1])));
        }
        
        // mapped image size
        long N[2] = {
            ceil(B[0] - A[0]) + 0.5,
            ceil(B[1] - A[1]) + 0.5
        };
        
        // pixels for mapped image
        const int nmap = N[0]*N[1];
        double* map = calloc(nmap, sizeof(double));
        if(!map)
            goto err_malloc;
        
        // set initally to null value if asked to
        if(use_null)
        {
            for(int j = 0; j < nmap; ++j)
                map[j] = NULVAL;
        }
        
        // map reference image onto image 1
        for(int k = 1, h = 0; k <= n[1]; ++k)
        for(int j = 1; j <= n[0]; ++j, ++h)
        {
            // pixel value
            const double f = (1./NSUB/NSUB)*det*pix[h];
            
            // subsampling
            for(long t = 0; t < NSUB; ++t)
            for(long s = 0; s < NSUB; ++s)
            {
                // position in reference image coordinate system
                const double X[2] = {
                    j - 0.5 + (s + 0.5)/NSUB - a[0],
                    k - 0.5 + (t + 0.5)/NSUB - a[1]
                };
                
                // position in image coordinate system
                const double Y[2] = {
                    a[2*i + 0] + T[0]*X[0] + T[1]*X[1] + 1 - A[0],
                    a[2*i + 1] + T[2]*X[0] + T[3]*X[1] + 1 - A[1]
                };
                
                // set pixel if within bounds
                if(Y[0] >= 0.5 && Y[0] < N[0] + 0.5 &&
                   Y[1] >= 0.5 && Y[1] < N[1] + 0.5)
                {
                    // convert to array indices
                    const int u = Y[0] - 0.5;
                    const int v = Y[1] - 0.5;
                    
                    // mapped pixel index
                    const int m = N[0]*v + u;
                    
                    // activate mapped pixel if necessary
                    if(map[m] == NULVAL)
                        map[m] = 0;
                    
                    // map pixel onto image
                    map[m] += f;
                }
            }
        }
        
        // create image extension for mapped image
        fits_create_img(fp, DOUBLE_IMG, 2, N, &err);
        
        // metadata for FITS
        if(i == 0)
        {
            // record file origin
            fits_write_key(fp, TSTRING, "ORIGIN", "immap",
                           "FITS file originator", &err);
            
            // record the date of FITS creation
            fits_write_date(fp, &err);
        }
        
        // set extension name
        snprintf(extname, sizeof(extname), "IMG%d", i);
        fits_write_key(fp, TSTRING, "EXTNAME", extname, "extension name",
                       &err);
        
        // write mapped image
        fits_write_pixnull(fp, TDOUBLE, fpix, nmap, map, nulv, &err);
        
        // done with mapped image
        free(map);
        
        // check for errors
        if(err)
            goto err_outfits;
    }
    
    // close output FITS
    fits_close_file(fp, &err);
    
    // check for errors
    if(err)
        goto err_outfits;
    
    // clean up
    free(pix);
    free(mat);
    if(x)
        free(x);
    free(a);
    
    
    // done
    return EXIT_SUCCESS;
    
    // errors
    
err_usage:
    fprintf(err ? stderr : stdout, "%s\n", USAGE);
    if(err && msg)
        fprintf(stderr, "\nerror: %s\n", msg);
    return err ? EXIT_FAILURE : EXIT_SUCCESS;
    
err_malloc:
    perror(msg);
    return EXIT_FAILURE;
    
err_matfile:
    fprintf(stderr, "%s: %s\n", matfile, msg);
    return EXIT_FAILURE;
    
err_ptsfile:
    fprintf(stderr, "%s: %s\n", ptsfile, msg);
    return EXIT_FAILURE;
    
err_ancfile:
    fprintf(stderr, "%s: %s\n", ancfile, msg);
    return EXIT_FAILURE;
    
err_imfits:
    msg = imfile;
    goto err_fitsio;
    
err_outfits:
    msg = outfile;
    goto err_fitsio;
    
err_fitsio:
    if(err)
    {
        char status[32];
        fits_get_errstatus(err, status);
        fprintf(stderr, "%s: %s\n", msg, status);
    }
    return EXIT_FAILURE;
}

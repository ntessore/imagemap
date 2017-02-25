#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "fitsio.h"

#include "input.h"

static const char* USAGE =
"usage: srcim [-0c] [-p PTSFILE] [-a ANCFILE] [-d PAD] \n"
"             OUTFITS INFILE IMFITS0 [IMFITS1 ...]";

// subsampling
#ifndef NSUB
#define NSUB 32
#endif

// value for null pixels
static const double NULVAL = DOUBLENULLVALUE;

int main(int argc, char* argv[])
{
    int err = 0;
    const char* msg = NULL;
    
    // input
    char* outfile;
    char* infile;
    char* ptsfile;
    char* ancfile;
    int ni;
    const char** images;
    int use_null, combine;
    double pad;
    
    // points
    int nx;
    double* x;
    double* a;
    
    // convergence ratio and shear table
    int nrow, ncol;
    double* fg;
    
    // FITS file
    fitsfile* infits;
    fitsfile* outfits;
    int na;
    long n[2], npix, fpix[2] = { 1, 1 };
    void* nulv;
    char extname[10];
    double* pix;
    
    // source image bounding box
    double A[2], B[2];
    
    // source image
    long N[2];
    long nsrc;
    double* src;
    double* wht;
    
    // default arguments
    outfile = infile = ptsfile = ancfile = NULL;
    ni = 0;
    images = NULL;
    use_null = combine = 0;
    pad = 0.;
    
    // parse arguments
    err = 0;
    for(int i = 1; i < argc && !err; ++i)
    {
        if(argv[i][0] == '-' && !(argv[i][1] == '-' && argv[i][2] == '\0'))
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
                // combine source images
                else if(*c == 'c')
                {
                    if(!combine)
                        combine = 1;
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
                // source image padding
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
            if(!outfile)
                outfile = argv[i];
            else if(!infile)
                infile = argv[i];
            else
            {
                images = realloc(images, (ni+1)*sizeof(double*));
                if(!images)
                    goto err_malloc;
                images[ni++] = argv[i];
            }
        }
    }
    
    // make sure input was given
    if(!outfile || !infile || !ni)
        err = 1;
    
    // combine mode needs points
    if(combine && !ptsfile)
    {
        msg = "combine mode (-c) needs points (-p)";
        err = 1;
    }
    
    // anchor points need combine mode
    if(ancfile && !combine)
    {
        msg = "anchor points (-a) can only be given in combine mode (-c)";
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
    
    // read points if given
    if(ptsfile)
    {
        // number of images from points file
        int nix;
        
        // read points from file
        read_points(ptsfile, &nix, &nx, &x);
        
        // if single image was given, use points file for image specification
        if(ni == 1 && nix > 1)
        {
            ni = nix;
            images = realloc(images, ni*sizeof(double*));
            if(!images)
                goto err_malloc;
            for(int i = 1; i < ni; ++i)
                images[i] = images[0];
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
        read_table(ancfile, &nrow, &ncol, &a);
        
        // check format
        if(ncol != 2)
        {
            msg = "anchor points file needs two columns (x, y)";
            goto err_ancfile;
        }
        
        // check number of images
        if(nrow != ni)
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
    
    // read input file
    read_table(infile, &nrow, &ncol, &fg);
    
    // make sure format is correct
    if(ncol != 3)
    {
        msg = "input file needs three columns (f, g1, g2)";
        goto err_infile;
    }
    
    // make sure number of images matches number of rows
    if(nrow != ni)
    {
        msg = "number of images must match number of rows";
        goto err_infile;
    }
    
    // create FITS file
    fits_create_file(&outfits, outfile, &err);
    
    // check for FITS creation errors
    if(err)
    {
        msg = outfile;
        goto err_fitsio;
    }
    
    // start with no source or weight image
    src = wht = NULL;
    
    // transform and write images one by one
    for(int i = 0; i < ni; ++i)
    {
        // convergence ratio and shear of image 0 and image 1
        const double f  = fg[3*0+0];
        const double g1 = fg[3*0+1];
        const double g2 = fg[3*0+2];
        const double F  = fg[3*i+0];
        const double G1 = fg[3*i+1];
        const double G2 = fg[3*i+2];
        
        // magnification matrix normalisation
        const double K = f/F/sqrt(fabs(1 - g1*g1 - g2*g2));
        
        // normalised magnification matrix for image i
        const double M[4] = { K*(1 - G1), K*(-G2), K*(-G2), K*(1 + G1) };
        
        // determinant of magnification matrix
        const double det = fabs(M[0]*M[3] - M[1]*M[2]);
        
        // read image if not skipped
        if(strcmp(images[i], "--") != 0)
        {
            // open input FITS file
            fits_open_image(&infits, images[i], READONLY, &err);
            fits_get_img_dim(infits, &na, &err);
            fits_get_img_size(infits, 2, n, &err);
            
            // check for errors
            if(err)
            {
                msg = images[i];
                goto err_fitsio;
            }
            
            // make sure the image is an image
            if(na != 2)
            {
                err = BAD_DIMEN;
                msg = images[i];
                goto err_fitsio;
            }
            
            // array for image pixels
            npix = n[0]*n[1];
            pix = malloc(npix*sizeof(double));
            if(!pix)
            {
                msg = images[i];
                goto err_malloc;
            }
            
            // read pixels into array and close file
            fits_read_pix(infits, TDOUBLE, fpix, npix, NULL, pix, NULL, &err);
            fits_close_file(infits, &err);
            
            // check for errors
            if(err)
            {
                msg = images[i];
                goto err_fitsio;
            }
        }
        else
        {
            // no image read
            npix = n[0] = n[1] = 0;
        }
        
        // create source plane if necessary
        if(!src)
        {
            // get bounding box of source plane
            if(nx > 0)
            {
                // get bounding box by mapping points
                for(int j = 0; j < nx; ++j)
                {
                    // get point in image plane, relative to anchor point
                    const double p[2] = {
                        x[5*nx*i + 5*j + 0] - a[2*i + 0],
                        x[5*nx*i + 5*j + 1] - a[2*i + 1]
                    };
                    
                    // transform point to source plane
                    const double q[2] = {
                        M[0]*p[0] + M[1]*p[1],
                        M[2]*p[0] + M[3]*p[1]
                    };
                    
                    // expand bounding box if necessary
                    if(j == 0 || M[0] > q[0])
                        A[0] = q[0];
                    if(j == 0 || A[1] > q[1])
                        A[1] = q[1];
                    if(j == 0 || B[0] < q[0])
                        B[0] = q[0];
                    if(j == 0 || B[1] < q[1])
                        B[1] = q[1];
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
            else if(n[0] > 0 && n[1] > 0)
            {
                // map image corners
                const double C[4][2] = {
                    { M[0]*1    + M[1]*1   , M[2]*1    + M[3]*1    },
                    { M[0]*n[0] + M[1]*1   , M[2]*n[0] + M[3]*1    },
                    { M[0]*1    + M[1]*n[1], M[2]*1    + M[3]*n[1] },
                    { M[0]*n[0] + M[1]*n[1], M[2]*n[0] + M[3]*n[1] }
                };
                
                // get source image bounding box
                A[0] = fmin(C[0][0], fmin(C[1][0], fmin(C[2][0], C[3][0])));
                A[1] = fmin(C[0][1], fmin(C[1][1], fmin(C[2][1], C[3][1])));
                B[0] = fmax(C[0][0], fmax(C[1][0], fmax(C[2][0], C[3][0])));
                B[1] = fmax(C[0][1], fmax(C[1][1], fmax(C[2][1], C[3][1])));
            }
            else
            {
                // no bounding box possible
                A[0] = B[0] = 0;
                A[1] = B[1] = 0;
            }
            
            // create source plane if there is a bounding box
            if(A[0] != B[0] && A[1] != B[1])
            {
                // get source image size
                N[0] = ceil(B[0] - A[0]) + 0.5;
                N[1] = ceil(B[1] - A[1]) + 0.5;
                
                // total number of source pixels
                nsrc = N[0]*N[1];
                
                // create array for source pixels and weights
                src = calloc(nsrc, sizeof(double));
                wht = calloc(nsrc, sizeof(double));
                if(!src || !wht)
                    goto err_malloc;
            }
            else
            {
                // no source plane
                N[0] = N[1] = 0;
            }
        }
        
        // map image pixels to source plane if possible
        if(npix > 0 && src)
        {
            for(int k = 1, h = 0; k <= n[1]; ++k)
            for(int j = 1; j <= n[0]; ++j, ++h)
            {
                // pixel weight and value
                const double w = (1./NSUB/NSUB)*det;
                const double f = w*pix[h];
                
                // subsampling
                for(long t = 0; t < NSUB; ++t)
                for(long s = 0; s < NSUB; ++s)
                {
                    // position in image coordinate system
                    const double dx = j - 0.5 + (s + 0.5)/NSUB - a[2*i + 0];
                    const double dy = k - 0.5 + (t + 0.5)/NSUB - a[2*i + 1];
                    
                    // position in source coordinate system
                    const double x = M[0]*dx + M[1]*dy - A[0];
                    const double y = M[2]*dx + M[3]*dy - A[1];
                    
                    // convert to array indices
                    const int u = x - 0.5;
                    const int v = y - 0.5;
                    
                    // map to source plane if within bounds
                    if(x >= 0.5 && x < N[0] + 0.5 &&
                       y >= 0.5 && y < N[1] + 0.5)
                    {
                        // index of source pixel
                        const int s = v*N[0]+u;
                        
                        // map pixel to source
                        src[s] += f;
                        wht[s] += w;
                    }
                }
            }
            
            // done with pixels
            free(pix);
        }
        
        // create FITS extension if not in combine mode
        if(src && !combine)
        {
            // create image extension for source in output file
            fits_create_img(outfits, DOUBLE_IMG, 2, N, &err);
            
            // record file origin
            fits_write_key(outfits, TSTRING, "ORIGIN", "srcim",
                           "FITS file originator", &err);
            
            // record the date of FITS creation
            fits_write_date(outfits, &err);
            
            // set extension name
            snprintf(extname, sizeof(extname), "SRC%d", i);
            fits_write_key(outfits, TSTRING, "EXTNAME", extname, 
                           "extension name", &err);
            
            // write source pixels if present and destroy source plane
            if(src)
            {
                // set unmapped source pixels to null value if asked to
                if(use_null)
                {
                    for(int j = 0; j < nsrc; ++j)
                        if(!wht[j])
                            src[j] = NULVAL;
                }
                
                // write source image
                fits_write_pixnull(outfits, TDOUBLE, fpix, nsrc, src, nulv,
                                   &err);
                
                // destroy source plane
                free(src);
                free(wht);
                src = wht = NULL;
            }
            
            // check for errors
            if(err)
            {
                msg = outfile;
                goto err_fitsio;
            }
        }
    }
    
    // write source plane in combine mode
    if(combine)
    {
        // normalise source pixels by total weight
        for(int i = 0; i < nsrc; ++i)
            if(wht[i] > 0)
                src[i] /= wht[i];
        
        // create image extension for source in output file
        fits_create_img(outfits, DOUBLE_IMG, 2, N, &err);
        
        // record file origin
        fits_write_key(outfits, TSTRING, "ORIGIN", "srcim",
                       "FITS file originator", &err);
        
        // record the date of FITS creation
        fits_write_date(outfits, &err);
        
        // set extension name
        fits_write_key(outfits, TSTRING, "EXTNAME", "SRC", "extension name",
                       &err);
        
        // write source image
        fits_write_pixnull(outfits, TDOUBLE, fpix, nsrc, src, nulv, &err);
        
        // create image extension for weights in output file
        fits_create_img(outfits, DOUBLE_IMG, 2, N, &err);
        
        // set extension name
        fits_write_key(outfits, TSTRING, "EXTNAME", "WHT", "extension name",
                       &err);
        
        // write weight image
        fits_write_pixnull(outfits, TDOUBLE, fpix, nsrc, wht, nulv, &err);
        
        // free source pixels and weights
        free(src);
        free(wht);
    }
    
    // close output FITS file
    fits_close_file(outfits, &err);
    
    // check for errors
    if(err)
    {
        msg = outfile;
        goto err_fitsio;
    }
    
    // clean up
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
    
err_infile:
    fprintf(stderr, "%s: %s\n", infile, msg);
    return EXIT_FAILURE;
    
err_ancfile:
    fprintf(stderr, "%s: %s\n", ancfile, msg);
    return EXIT_FAILURE;
    
err_malloc:
    perror(msg);
    return EXIT_FAILURE;
    
err_fitsio:
    if(err)
    {
        char status[32];
        fits_get_errstatus(err, status);
        fprintf(stderr, "%s: %s\n", msg, status);
    }
    return EXIT_FAILURE;
}

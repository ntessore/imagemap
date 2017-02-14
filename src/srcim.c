#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "fitsio.h"

#include "input.h"

// subsampling
#ifndef NSUB
#define NSUB 32
#endif

// report FITS error and exit
void fits_error(const char* file, int status)
{
    char message[32];
    fits_get_errstatus(status, message);
    fprintf(stderr, "%s: %s\n", file, message);
    exit(1);
}

int main(int argc, char* argv[])
{
    int err;
    
    // input
    char* o;
    char* f;
    int nimg;
    const char** img;
    
    // convergence ratio and shear table
    int nrow, ncol;
    double* fg;
    
    // FITS file
    fitsfile* infits;
    fitsfile* outfits;
    int na;
    long n[2], npix, fpix[2] = { 1, 1 };
    char extname[10];
    double* pix;
    
    // source image bounding box
    double C[4][2], X[2], Y[2];
    
    // source image
    int NA;
    long N[2];
    long nsrc;
    double* src;
    
    // default arguments
    o = f = NULL;
    nimg = 0;
    
    // allocate space for max. number of possible arguments
    img = malloc(argc*sizeof(char*));
    if(!img)
    {
        perror(NULL);
        return EXIT_FAILURE;
    }
    
    // parse arguments
    err = 0;
    for(int i = 1; i < argc && !err; ++i)
    {
        if(argv[i][0] == '-' && !(argv[i][1] == '-' && argv[i][2] == '\0'))
        {
            // flags
            for(char* a = &argv[i][1]; *a; ++a)
            {
                // unknown flag
                {
                    err = 1;
                }
            }
        }
        else
        {
            // positional arguments
            if(!o)
                o = argv[i];
            else if(!f)
                f = argv[i];
            else
                img[nimg++] = argv[i];
        }
    }
    
    // make sure input was given
    if(!o || !f || !nimg)
        err = 1;
    
    // check for input errors
    if(err)
    {
        fprintf(stderr, "usage: srcim OUTFITS INFILE IMFITS0 [IMFITS1 ...]\n");
        return EXIT_FAILURE;
    }
    
    // read input file
    read_table(f, &nrow, &ncol, &fg);
    
    // make sure format is correct
    if(ncol != 3)
    {
        fprintf(stderr, "%s: input file needs three columns (f, g1, g2)\n", f);
        return EXIT_FAILURE;
    }
    
    // make sure number of images matches number of rows
    if(nrow != nimg)
    {
        fprintf(stderr, "%s: number of images must match number of rows\n", f);
        return EXIT_FAILURE;
    }
    
    // create FITS file
    fits_create_file(&outfits, o, &err);
    
    // check for FITS creation errors
    if(err)
        fits_error(o, err);
    
    // transform and write images one by one
    for(int i = 0; i < nimg; ++i)
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
        const double A[4] = { K*(1 - G1), K*(-G2), K*(-G2), K*(1 + G1) };
        
        // determinant of magnification matrix
        const double det = fabs(A[0]*A[3] - A[1]*A[2]);
        
        // check whether to skip image
        if(strcmp(img[i], "--") == 0)
        {
            // empty image extension for this source image
            NA = 0;
            N[0] = 0;
            N[1] = 0;
        }
        else
        {
            // open input FITS file
            fits_open_image(&infits, img[i], READONLY, &err);
            fits_get_img_dim(infits, &na, &err);
            fits_get_img_size(infits, 2, n, &err);
            
            // check for errors
            if(err)
                fits_error(img[i], err);
            
            // make sure the image is an image
            if(na != 2)
                fits_error(img[i], BAD_DIMEN);
            
            // array for image pixels
            npix = n[0]*n[1];
            pix = malloc(npix*sizeof(double));
            if(!pix)
            {
                perror(img[i]);
                return EXIT_FAILURE;
            }
            
            // read pixels into array and close file
            fits_read_pix(infits, TDOUBLE, fpix, npix, NULL, pix, NULL, &err);
            fits_close_file(infits, &err);
            
            // check for errors
            if(err)
                fits_error(img[i], err);
            
            // map image corners
            C[0][0] = 0;
            C[0][1] = 0;
            C[1][0] = A[0]*n[0];
            C[1][1] = A[2]*n[0];
            C[2][0] = A[1]*n[1];
            C[2][1] = A[3]*n[1];
            C[3][0] = A[0]*n[0] + A[1]*n[1];
            C[3][1] = A[2]*n[0] + A[3]*n[1];
            
            // get source image bounding box
            X[0] = fmin(C[0][0], fmin(C[1][0], fmin(C[2][0], C[3][0])));
            X[1] = fmin(C[0][1], fmin(C[1][1], fmin(C[2][1], C[3][1])));
            Y[0] = fmax(C[0][0], fmax(C[1][0], fmax(C[2][0], C[3][0])));
            Y[1] = fmax(C[0][1], fmax(C[1][1], fmax(C[2][1], C[3][1])));
            
            // two dimensions for image
            NA = 2;
            
            // get source image size
            N[0] = ceil(Y[0] - X[0]) + 0.5;
            N[1] = ceil(Y[1] - X[1]) + 0.5;
            
            // array for source pixels
            nsrc = N[0]*N[1];
            src = calloc(nsrc, sizeof(double));
            if(!src)
            {
                perror(img[i]);
                return EXIT_FAILURE;
            }
            
            // map image to source plane
            for(int j = 0, k = 0; j < n[1]; ++j)
            for(int i = 0; i < n[0]; ++i, ++k)
            {
                // pixel value
                const double f = (1./NSUB/NSUB)*det*pix[k];
                
                // subsampling
                for(long t = 0; t < NSUB; ++t)
                for(long s = 0; s < NSUB; ++s)
                {
                    // position in image coordinate system
                    const double dx = i - 0.5 + (s + 0.5)/NSUB;
                    const double dy = j - 0.5 + (t + 0.5)/NSUB;
                    
                    // position in source coordinate system
                    const double x = A[0]*dx + A[1]*dy - X[0];
                    const double y = A[2]*dx + A[3]*dy - X[1];
                    
                    // convert to array indices
                    const int u = x - 0.5;
                    const int v = y - 0.5;
                    
                    // ignore out of bounds pixels
                    if(x < 0.5 || x >= N[0] + 0.5 ||
                       y < 0.5 || y >= N[1] + 0.5)
                        continue;
                    
                    // map pixel to source
                    src[v*N[0]+u] += f;
                }
            }
        }
        
        // create image extension for source in output file
        fits_create_img(outfits, DOUBLE_IMG, NA, N, &err);
        
        // metadata for primary HDU
        if(i == 0)
        {
            // record file origin
            fits_write_key(outfits, TSTRING, "ORIGIN", "srcim",
                           "FITS file originator", &err);
            
            // record the date of FITS creation
            fits_write_date(outfits, &err);
        }
        
        // set extension name
        snprintf(extname, sizeof(extname), "SRC%d", i);
        fits_write_key(outfits, TSTRING,
                       "EXTNAME", extname, "extension name", &err);
        
        // process source pixels if not skipped
        if(NA != 0)
        {
            // write source image
            fits_write_pix(outfits, TDOUBLE, fpix, nsrc, src, &err);
            
            // done with pixels
            free(pix);
            free(src);
        }
        
        // check for errors
        if(err)
            fits_error(o, err);
    }
    
    // close output FITS file
    fits_close_file(outfits, &err);
    
    // check for errors
    if(err)
        fits_error(o, err);
    
    return EXIT_SUCCESS;
}

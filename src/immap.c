#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fitsio.h"

#include "input.h"

const char* USAGE =
"usage: immap MATFILE IMFITS OUTFITS";

// subsampling
#ifndef NSUB
#define NSUB 32
#endif

int main(int argc, char* argv[])
{
    // error flag and message
    int err = 0;
    const char* msg = NULL;
    
    // input
    char* matfile;
    char* imfile;
    char* outfile;
    
    // matrix table
    int nrow, ncol;
    double* T;
    
    // FITS file
    fitsfile* fp;
    int na;
    long n[2], npix, fpix[2] = { 1, 1 };
    double* pix;
    
    // default arguments
    matfile = imfile = outfile = NULL;
    
    // parse arguments
    err = 0;
    for(int i = 1; i < argc && !err; ++i)
    {
        if(argv[i][0] == '-')
        {
            // flags
            for(char* c = &argv[i][1]; *c; ++c)
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
    
    // check for input errors
    if(err)
        goto err_usage;
    
    // read input file
    read_table(matfile, &nrow, &ncol, &T);
    
    // make sure format is correct
    if(ncol != 4)
    {
        msg = "file needs four columns (matrix entries)";
        goto err_matfile;
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
    
    // create image extension for reference image
    fits_create_img(fp, DOUBLE_IMG, 2, n, &err);
    
    // record file origin
    fits_write_key(fp, TSTRING, "ORIGIN", "immap", "FITS file originator",
                   &err);
    
    // record the date of FITS creation
    fits_write_date(fp, &err);
    
    // write reference image
    fits_write_pix(fp, TDOUBLE, fpix, npix, pix, &err);
    
    // check for errors
    if(err)
        goto err_outfits;
    
    // transform and write images
    for(int i = 0; !err && i < nrow; ++i)
    {
        // matrix coefficients
        const double A = T[4*i+0], B = T[4*i+1], C = T[4*i+2], D = T[4*i+3];
        
        // determinant of matrix
        const double det = fabs(A*D - B*C);
        
        // map image corners
        const double c[4][2] = {
            {               0,               0 },
            {          A*n[0],          C*n[0] },
            {          B*n[1],          D*n[1] },
            { A*n[0] + B*n[1], C*n[0] + D*n[1] }
        };
        
        // get mapped image bounding box
        const double X[2] = {
            fmin(c[0][0], fmin(c[1][0], fmin(c[2][0], c[3][0]))),
            fmin(c[0][1], fmin(c[1][1], fmin(c[2][1], c[3][1]))),
        };
        const double Y[2] = {
            fmax(c[0][0], fmax(c[1][0], fmax(c[2][0], c[3][0]))),
            fmax(c[0][1], fmax(c[1][1], fmax(c[2][1], c[3][1]))),
        };
        
        // mapped image size
        long N[2] = { ceil(Y[0]-X[0]) + 0.5, ceil(Y[1]-X[1]) + 0.5 };
        
        // pixels for mapped image
        const int nmap = N[0]*N[1];
        double* map = calloc(nmap, sizeof(double));
        if(!map)
            goto err_malloc;
        
        // map reference image onto image 1
        for(int j = 0, k = 0; j < n[1]; ++j)
        for(int i = 0; i < n[0]; ++i, ++k)
        {
            // pixel value
            const double f = (1./NSUB/NSUB)*det*pix[k];
            
            // subsampling
            for(long t = 0; t < NSUB; ++t)
            for(long s = 0; s < NSUB; ++s)
            {
                // position in reference image coordinate system
                const double dx = i - 0.5 + (s + 0.5)/NSUB;
                const double dy = j - 0.5 + (t + 0.5)/NSUB;
                
                // position in image coordinate system
                const double x = A*dx + B*dy - X[0];
                const double y = C*dx + D*dy - X[1];
                
                // convert to array indices
                const int u = x - 0.5;
                const int v = y - 0.5;
                
                // ignore out of bounds pixels
                if(x < 0.5 || x >= N[0] + 0.5 || y < 0.5 || y >= N[1] + 0.5)
                    continue;
                
                // map pixel onto image
                map[v*N[0]+u] += f;
            }
        }
        
        // create image extension for mapped image
        fits_create_img(fp, DOUBLE_IMG, 2, N, &err);
        
        // write mapped image
        fits_write_pix(fp, TDOUBLE, fpix, nmap, map, &err);
        
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
    free(T);
    
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

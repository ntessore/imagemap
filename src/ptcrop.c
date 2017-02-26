#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fitsio.h"

#include "input.h"

const char* USAGE =
"usage: ptcrop [-0] [-d PAD] PTSFILE IMFITS OUTFITS";

// value for null pixels
static const double NULVAL = DOUBLENULLVALUE;

int main(int argc, char* argv[])
{
    // error flag and message
    int err = 0;
    const char* msg = NULL;
    
    // input
    char* ptsfile;
    char* imfile;
    char* outfile;
    int use_null;
    double pad;
    
    // points
    int ni, nx;
    double* x;
    
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
    ptsfile = imfile = outfile = NULL;
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
                // image padding
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
            if(!ptsfile)
                ptsfile = argv[i];
            else if(!imfile)
                imfile = argv[i];
            else if(!outfile)
                outfile = argv[i];
            else
                err = 1;
        }
    }
    
    // make sure input files were given
    if(!ptsfile || !imfile || !outfile)
        err = 1;
    
    // check for input errors
    if(err)
        goto err_usage;
    
    // set default padding to 50%
    if(!pad)
        pad = -0.5;
    
    // set null value if not writing zeros
    if(use_null)
        nulv = (void*)&NULVAL;
    else
        nulv = NULL;
    
    // read points from file
    read_points(ptsfile, &ni, &nx, &x);
    
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
    
    // crop images
    for(int i = 0; i < ni; ++i)
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
        
        // cropped image size
        long N[2] = {
            ceil(B[0] - A[0]) + 0.5,
            ceil(B[1] - A[1]) + 0.5
        };
        
        // pixels for cropped image
        const int ncrp = N[0]*N[1];
        double* crp = calloc(ncrp, sizeof(double));
        if(!crp)
            goto err_malloc;
        
        // set initally to null value if asked to
        if(use_null)
        {
            for(int j = 0; j < ncrp; ++j)
                crp[j] = NULVAL;
        }
        
        // map reference image onto image 1
        for(int k = 1, h = 0; k <= N[1]; ++k)
        for(int j = 1; j <= N[0]; ++j, ++h)
        {
            // position in original image
            const double x = A[0] - 1 + j;
            const double y = A[1] - 1 + k;
            
            // set pixel if within bounds
            if(x >= 0.5 && x < n[0] + 0.5 && y >= 0.5 && y < n[1] + 0.5)
            {
                // array indices
                const int u = x - 0.5;
                const int v = y - 0.5;
                const int w = n[0]*v + u;
                
                // copy pixel
                crp[h] = pix[w];
            }
        }
        
        // create image extension for cropped image
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
        
        // write cropped image
        fits_write_pixnull(fp, TDOUBLE, fpix, ncrp, crp, nulv, &err);
        
        // done with cropped image
        free(crp);
        
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
    if(x)
        free(x);
    
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

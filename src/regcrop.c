#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "fitsio.h"
#include "regions.h"

void error(const char* msg, int status)
{
    if(!msg)
        msg = strerror(errno);
    
    if(status)
    {
        // print FITS file and status error
        char buf[32] = {0};
        fits_get_errstatus(status, buf);
        fprintf(stderr, "%s: error: %s\n", msg, buf);
    }
    else
    {
        // print the message
        fprintf(stderr, "error: %s\n", msg);
    }
    
    exit(1);
}

int main(int argc, char* argv[])
{
    // input and output filenames
    const char* imgname = NULL;
    const char* regname = NULL;
    const char* outname = NULL;
    
    // FITS file reading and writing
    int status = 0;
    fitsfile* fptr;
    long naxis[2] = { 0, 0 };
    long fpixel[2] = { 1, 1 };
    char* fitshdr;
    int nkeys;
    int n1, n2;
    double* img;
    double* pix;
    
    // region data
    char* regions;
    int nmask, nreg;
    RegionsMask mask;
    int* sizes;
    
    // check for correct number of arguments
    if(argc != 4)
    {
        fprintf(stderr, "usage: regcrop FITSFILE REGIONFILE OUTFILE\n");
        exit(1);
    }
    
    // collect arguments
    imgname = argv[1];
    regname = argv[2];
    outname = argv[3];
    
    // open FITS image for reading and get size
    fits_open_image(&fptr, imgname, READONLY, &status);
    fits_get_img_size(fptr, 2, naxis, &status);
    
    // size of input
    n1 = naxis[0];
    n2 = naxis[1];
    
    // create array for input pixels
    img = malloc(n1*n2*sizeof(double));
    if(!img)
        error(NULL, 0);
    
    // read input pixels into array, read input header card, close FITS
    fits_read_pix(fptr, TDOUBLE, fpixel, n1*n2, NULL, img, NULL, &status);
    fits_convert_hdr2str(fptr, 0, NULL, 0, &fitshdr, &nkeys, &status);
    fits_close_file(fptr, &status);
    
    // check FITS operations
    if(status)
        error(imgname, status);
    
    // read regions from file
    regions = malloc(strlen(regname) + 2);
    if(!regions)
        error(NULL, 0);
    sprintf(regions, "@%s", regname);
    
    // open region definition
    Regions reg = OpenRegions(fitshdr, regions, NULL);
    
    // filter regions in image
    nmask = FilterRegions(reg, 1, n1, 1, n2, 1, &mask, &nreg);
    
    // initialise array of region sizes
    sizes = calloc(nreg*4, sizeof(int));
    if(!sizes)
        error(NULL, 0);
    
    // first pass: find region sizes
    for(int i = 0; i < nmask; ++i)
    {
        // size of current region
        int* s = sizes + 4*(mask[i].region-1);
        
        // get min and max extent
        if(!s[0] || mask[i].xstart < s[0])
            s[0] = mask[i].xstart;
        if(!s[1] || mask[i].xstop > s[1])
            s[1] = mask[i].xstop;
        if(!s[2] || mask[i].y < s[2])
            s[2] = mask[i].y;
        if(!s[3] || mask[i].y > s[3])
            s[3] = mask[i].y;
    }
    
    // create FITS for regions
    fits_create_file(&fptr, outname, &status);
    if(status)
        error(outname, status);
    
    // write primary HDU
    fits_create_img(fptr, SHORT_IMG, 0, NULL, &status);
    fits_write_key(fptr, TSTRING, "ORIGIN", "regcrop", "FITS file originator", &status);
    fits_write_date(fptr, &status);
    if(status)
        error(outname, status);
    
    // second pass: copy and write regions
    for(int i = 0; i < nreg; ++i)
    {
        // get region size
        int x0 = sizes[4*i+0];
        int x1 = sizes[4*i+1];
        int y0 = sizes[4*i+2];
        int y1 = sizes[4*i+3];
        
        // pixels in region
        naxis[0] = x1 - x0 + 1;
        naxis[1] = y1 - y0 + 1;
        
        // create array for region pixels
        pix = calloc(naxis[0]*naxis[1], sizeof(double));
        if(!pix)
            error(NULL, 0);
        
        // go through regions and copy pixels
        for(int j = 0; j < nmask; ++j)
        {
            // skip if not this region
            if(mask[j].region != i+1)
                continue;
            
            // copy pixels
            for(int x = mask[j].xstart, y = mask[j].y; x <= mask[j].xstop; ++x)
                pix[(y-y0)*naxis[0]+x-x0] = img[(y-1)*n1+x-1];
        }
        
        // create image extension for region
        fits_create_img(fptr, DOUBLE_IMG, 2, naxis, &status);
        if(status)
            error(outname, status);
        
        // write region pixels
        fits_write_pix(fptr, TDOUBLE, fpixel, naxis[0]*naxis[1], pix, &status);
        if(status)
            error(outname, status);
        
        // done with pixels
        free(pix);
    }
    
    // close region FITS file
    fits_close_file(fptr, &status);
    if(status)
        error(outname, status);
    
    // done
    CloseRegions(reg);
    free(mask);
    free(sizes);
    free(regions);
    free(fitshdr);
    free(img);
    
    return EXIT_SUCCESS;
}

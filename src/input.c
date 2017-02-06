#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "input.h"

// maximum length of input line
#ifndef LINELEN
#define LINELEN 1024
#endif

// newline characters
const char* TOK_NL = "\r\n";

// whitespace charaters
const char* TOK_WS = " \t";

// delimiter characters
const char* TOK_DL = ";";

void read_points(const char* filename, int* nimg, int* npts, double** pts)
{
    // file handle
    FILE* fp;
    
    // buffer for lines
    char linebuf[LINELEN];
    
    // current part of line
    size_t begin, len;
    
    // number of images
    int ni = 0;
    
    // number of points
    int np = 0;
    
    // points array
    double* p;
    
    // number of columns in line
    int nc;
    
    // open fitfile
    fp = fopen(filename, "r");
    if(!fp)
    {
        perror(filename);
        exit(1);
    }
    
    // points array must be NULL for first call to realloc
    p = NULL;
    
    // read fitfile line by line
    for(int line = 1; fgets(linebuf, sizeof(linebuf), fp); ++line)
    {
        // get whole line
        len = strcspn(linebuf, TOK_NL);
        
        // make sure line was read completely
        if(!linebuf[len] && !feof(fp))
        {
            fprintf(stderr, "%s: line %d: line too long\n", filename, line);
            exit(1);
        }
        
        // strip newline
        linebuf[len] = 0;
        
        // skip initial whitespace
        begin = strspn(linebuf, TOK_WS);
        
        // skip commented or empty lines
        if(linebuf[begin] == '#' || linebuf[begin] == '\0')
            continue;
        
        // find number of columns
        for(nc = 0; begin <= len; ++nc)
            begin += strcspn(linebuf + begin, TOK_DL) + 1;
        
        // make sure minimum number of columns is given
        if(nc < 3)
        {
            fprintf(stderr, "%s: line %d: expected format: point 1; point 2; "
                    "point 3 [ ; point 4 ... ]\n", filename, line);
            exit(1);
        }
        
        // store number of points if this is the first image
        if(!np)
            np = nc;
        
        // number of points must match previous rows
        if(nc != np)
        {
            fprintf(stderr, "%s: line %d: number of points does not match "
                    "previous images\n", filename, line);
            exit(1);
        }
        
        // expand points and files array for image
        p = realloc(p, 5*(ni+1)*np*sizeof(double));
        if(!p)
        {
            perror(NULL);
            exit(1);
        }
        
        // back to beginning of line
        begin = 0;
        
        // now process columns
        for(int i = 0; i < nc; ++i)
        {
            // point data
            double x, y, dx, dy, rho;
            
            // skip initial whitespace
            begin += strspn(linebuf + begin, TOK_WS);
            
            // find delimiter to end column
            len = strcspn(linebuf + begin, TOK_DL);
            
            // terminate column string
            linebuf[begin + len] = 0;
            
            // parse point
            switch(sscanf(linebuf + begin, "%lf,%lf,%lf,%lf,%lf", &x, &y, &dx, 
                          &dy, &rho))
            {
                case 2:
                    // no uncertainty given: assume 1px and uncorrelated
                    dx = 1;
                    dy = 1;
                    rho = 0;
                    break;
                    
                case 3:
                    // only dx given, assume dy same and uncorrelated
                    dy = dx;
                    rho = 0;
                    break;
                    
                case 4:
                    // assume uncorrelated
                    rho = 0;
                    break;
                    
                case 5:
                    // all is well
                    break;
                    
                default:
                    // invalid format
                    fprintf(stderr, "%s: line %d: column %d: point must be in "
                            "the form \"x, y\"\n", filename, line, i+1);
                    exit(1);
            }
            
            // store point
            p[5*(ni*np+i)+0] = x;
            p[5*(ni*np+i)+1] = y;
            p[5*(ni*np+i)+2] = dx;
            p[5*(ni*np+i)+3] = dy;
            p[5*(ni*np+i)+4] = rho;
            
            // skip past column
            begin += len + 1;
        }
        
        // new image stored
        ni += 1;
    }
    
    // close fitfile
    fclose(fp);
    
    // store number of images and points
    *nimg = ni;
    *npts = np;
    *pts = p;
}

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "input.h"

// maximum length of input line
#ifndef LINELEN
#define LINELEN 1024
#endif

// newline characters
const char* TOK_NL = "\r\n";

void read_points(const char* filename, int* nimg, int* npts, double** pts)
{
    // file handle
    FILE* fp;
    
    // current line in file
    int line;
    
    // buffer for lines
    char linebuf[LINELEN];
    
    // pointer to beginning of line
    char* b;
    
    // length of line
    size_t len;
    
    // number of images
    int ni = 0;
    
    // number of points
    int np = 0;
    
    // points array, must be NULL for first call to realloc
    double* p = NULL;
    
    // group size
    int ng = 0;
    
    // open points file
    fp = fopen(filename, "r");
    if(!fp)
    {
        perror(filename);
        exit(1);
    }
    
    // read file line by line
    for(line = 1; fgets(linebuf, sizeof(linebuf), fp); ++line)
    {
        // point data
        double x, y, dx, dy, rho;
        
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
        b = linebuf;
        while(*b && isspace(*b))
            b += 1;
        
        // skip commented lines
        if(*b == '#')
            continue;
        
        // empty lines end groups of points
        if(*b == '\0')
        {
            // ignore blank lines between images
            if(np == ni*ng)
                continue;
            
            // first group of points defines number of points per image
            if(!ng)
                ng = np;
            
            // new image stored
            ni += 1;
            
            // number of points must match previous images
            if(np != ni*ng)
                goto group_error;
            
            // done with line
            continue;
        }
        
        // parse point, with fallthrough
        switch(sscanf(b, "%lf%lf%lf%lf%lf", &x, &y, &dx, &dy, &rho))
        {
            case 2:
                // no uncertainty given: assume 1px and uncorrelated
                dx = 1;
                
            case 3:
                // only dx given, assume dy same and uncorrelated
                dy = dx;
                
            case 4:
                // assume uncorrelated
                rho = 0;
                
            case 5:
                // all is well
                break;
                
            default:
                // invalid format
                fprintf(stderr, "%s: line %d: point must be in the form "
                        "x y [dx [dy [rho]]]\n", filename, line);
                exit(1);
        }
        
        // expand points array
        p = realloc(p, 5*(np+1)*sizeof(double));
        if(!p)
        {
            perror(NULL);
            exit(1);
        }
        
        // store point
        p[5*np+0] = x;
        p[5*np+1] = y;
        p[5*np+2] = dx;
        p[5*np+3] = dy;
        p[5*np+4] = rho;
        
        // new point stored
        np += 1;
    }
    
    // close file
    fclose(fp);
    
    // check if there was a final, unterminated group in the file
    if(np != ni*ng)
    {
        if(np == (ni+1)*ng)
            ni += 1;
        else
            goto group_error;
    }
    
    // store number of images and points per image
    *nimg = ni;
    *npts = ng;
    *pts = p;
    
    // done
    return;
    
group_error:
    fprintf(stderr, "%s: line %d: number of points does not match "
            "previous images\n", filename, line);
    exit(1);
}

void read_table(const char* filename, int* nrow, int* ncol, double** tab)
{
    // file handle
    FILE* fp;
    
    // buffer for lines
    char linebuf[LINELEN];
    
    // current cell
    char* beg;
    char* end;
    
    // number of rows
    int nr = 0;
    
    // number of columns
    int nc = 0;
    
    // table array
    double* t;
    
    // open table file
    fp = fopen(filename, "r");
    if(!fp)
    {
        perror(filename);
        exit(1);
    }
    
    // table array must be NULL for first call to realloc
    t = NULL;
    
    // read file line by line
    for(int line = 1; fgets(linebuf, sizeof(linebuf), fp); ++line)
    {
        // get whole line
        end = linebuf + strcspn(linebuf, TOK_NL);
        
        // make sure line was read completely
        if(*end =='\0' && !feof(fp))
        {
            fprintf(stderr, "%s: line %d: line too long\n", filename, line);
            exit(1);
        }
        
        // strip newline
        *end = '\0';
        
        // skip initial whitespace
        beg = linebuf;
        while(*beg && isspace(*beg))
            ++beg;
        
        // skip commented or empty lines
        if(*beg == '#' || *beg == '\0')
            continue;
        
        // find number of columns if this is the first row
        if(!nc)
        {
            char* b = beg;
            
            for(;; ++nc)
            {
                // skip initial whitespace
                while(*b && isspace(*b))
                    ++b;
                
                // check if line ended
                if(*b == '\0')
                    break;
                
                // parse number
                strtod(b, &end);
                
                // check parsing
                if(*end && !isspace(*end))
                {
                    fprintf(stderr, "%s: line %d: column %d: parse error",
                            filename, line, nc+1);
                    exit(1);
                }
                
                // next column
                b = end;
            }
        }
        
        // expand table array for row
        t = realloc(t, (nr+1)*nc*sizeof(double));
        if(!t)
        {
            perror(NULL);
            exit(1);
        }
        
        // read columns
        for(int i = 0; i < nc; ++i)
        {
            // skip initial whitespace
            while(*beg && isspace(*beg))
                ++beg;
            
            // check if there is line left
            if(*beg == '\0')
            {
                fprintf(stderr, "%s: line %d: too few values", filename, line);
                exit(1);
            }
            
            // parse number
            t[nr*nc+i] = strtod(beg, &end);
            
            // check parsing
            if(*end && !isspace(*end))
            {
                fprintf(stderr, "%s: line %d: column %d: parse error",
                        filename, line, i+1);
                exit(1);
            }
            
            // next column
            beg = end;
        }
        
        // skip whitespace
        while(*beg && isspace(*beg))
            ++beg;
        
        // check if there is line left
        if(*beg != '\0')
        {
            fprintf(stderr, "%s: line %d: too many values", filename, line);
            exit(1);
        }
        
        // new row stored
        nr += 1;
    }
    
    // close file
    fclose(fp);
    
    // store table
    *nrow = nr;
    *ncol = nc;
    *tab = t;
}

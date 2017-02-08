#include <stdlib.h>
#include <stdio.h>

#include "input.h"

int main(int argc, char* argv[])
{
    int err;
    
    // input
    char* f;
    char* o;
    int v;
    
    // default arguments
    f = o = NULL;
    v = 0;
    
    // rows and columns of input file
    int nrow, ncol;
    double* tab;
    
    // output file pointer
    FILE* fp;
    
    // reference shear
    double g1, g2;
    
    // parse arguments
    err = 0;
    for(int i = 1; i < argc && !err; ++i)
    {
        if(argv[i][0] == '-')
        {
            // flags
            for(char* a = &argv[i][1]; *a; ++a)
            {
                // decrease verbosity
                if(*a == 'q')
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
        fprintf(stderr, "usage: mat2lens [-q] [-o OUTFILE] FILE\n");
        return EXIT_FAILURE;
    }
    
    // read input file
    read_table(f, &nrow, &ncol, &tab);
    
    // make sure format is correct
    if(ncol != 4)
    {
        fprintf(stderr, "%s: file needs four columns (matrix entries)\n", f);
        return EXIT_FAILURE;
    }
    
    // check that enough matrices are given
    if(nrow < 3)
    {
        fprintf(stderr, "%s: needs at least two matrices\n", f);
        return EXIT_FAILURE;
    }
    
    // open output file if given
    if(o)
    {
        fp = fopen(o, "w");
        if(!fp)
        {
            perror(o);
            return EXIT_FAILURE;
        }
    }
    else
    {
        // no output file
        fp = NULL;
    }
    
    // transform matrices to abcd coefficients
    for(int i = 0; i < nrow; ++i)
    {
        const double a = tab[4*i+0] - tab[4*i+3];
        const double b = tab[4*i+2] + tab[4*i+1];
        const double c = tab[4*i+2] - tab[4*i+1];
        const double d = tab[4*i+0] + tab[4*i+3];
        
        tab[4*i+0] = a;
        tab[4*i+1] = b;
        tab[4*i+2] = c;
        tab[4*i+3] = d;
    }
    
    // compute reference shear from matrix 0 and 1
    {
        const double a = tab[4*0+0];
        const double b = tab[4*0+1];
        const double c = tab[4*0+2];
        const double A = tab[4*1+0];
        const double B = tab[4*1+1];
        const double C = tab[4*1+2];
        
        // compute g for reference image 0
        g1 = (a*C - A*c)/(b*A - B*a);
        g2 = (b*C - B*c)/(b*A - B*a);
    }
    
    // output convergence ratio and shear of reference image
    if(v >= 0)
        printf("% 18.8f % 18.8f % 18.8f\n", 1.0, g1, g2);
    
    // compute f, g1, g2 for rest of images
    for(int i = 0; i < nrow; ++i)
    {
        // a,b,c,d coefficients for image
        const double a = tab[4*i+0];
        const double b = tab[4*i+1];
        const double c = tab[4*i+2];
        const double d = tab[4*i+3];
        
        // numerator and denominator for f and g
        const double J = 0.5*(c*c + d*d - a*a - b*b);
        const double P = d*g1 - c*g2 + a;
        const double Q = c*g1 + d*g2 + b;
        const double R = a*g1 + b*g2 + d;
        
        // f, g1, g2 for image
        const double F  = J/R;
        const double G1 = P/R;
        const double G2 = Q/R;
        
        // write convergence ratio and shear
        if(v >= 0)
            printf("% 18.8f % 18.8f % 18.8f\n", F, G1, G2);
        if(fp)
            fprintf(fp, "% 18.8f % 18.8f % 18.8f\n", F, G1, G2);
    }
    
    // close output file if open
    if(fp)
        fclose(fp);
    
    // done
    free(tab);
    
    return EXIT_SUCCESS;
}

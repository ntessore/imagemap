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
        fprintf(stderr, "usage: lens2mat [-q] [-o OUTFILE] FILE\n");
        return EXIT_FAILURE;
    }
    
    // read input file
    read_table(f, &nrow, &ncol, &tab);
    
    // make sure format is correct
    if(ncol != 3)
    {
        fprintf(stderr, "%s: file needs three columns (f, g1, g2)\n", f);
        return EXIT_FAILURE;
    }
    
    // check that enough images are given
    if(nrow < 2)
    {
        fprintf(stderr, "%s: needs at least two images\n", f);
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
    
    // compute transformation matrices
    for(int i = 1; i < nrow; ++i)
    {
        // convergence ratios and shears of images 0 and i
        const double f  = tab[0];
        const double g1 = tab[1];
        const double g2 = tab[2];
        const double F  = tab[3*i+0];
        const double G1 = tab[3*i+1];
        const double G2 = tab[3*i+2];
        
        // matrix coefficients
        const double
        A = (F/f)*((1 - g1)*(1 + G1) - g2*G2)/(1 - G1*G1 - G2*G2),
        B = (F/f)*((1 + g1)*G2 - (1 + G1)*g2)/(1 - G1*G1 - G2*G2),
        C = (F/f)*((1 - g1)*G2 - (1 - G1)*g2)/(1 - G1*G1 - G2*G2),
        D = (F/f)*((1 + g1)*(1 - G1) - g2*G2)/(1 - G1*G1 - G2*G2);
        
        // write matrix coefficients
        if(v >= 0)
            printf("% 18.8f % 18.8f % 18.8f % 18.8f\n", A, B, C, D);
        if(fp)
            fprintf(fp, "% 18.8f % 18.8f % 18.8f % 18.8f\n", A, B, C, D);
    }
    
    // close output file if open
    if(fp)
        fclose(fp);
    
    // done
    free(tab);
    
    return EXIT_SUCCESS;
}

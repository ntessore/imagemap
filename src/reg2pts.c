#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef LINELEN
#define LINELEN 1024
#endif

const char* TOK_NL = "\r\n";
const char* TOK_WS = " \t";
const char* TOK_DL = " (";

const double DEGREE = 0.017453292519943295769;

const double SCALE = 0.40853898265363493526;

int main(int argc, char* argv[])
{
    int err;
    
    int filc, f, n, m;
    const char** filv;
    
    FILE* fp;
    int line;
    char linebuf[LINELEN];
    size_t begin, len, tok;
    
    unsigned int code;
    
    // default arguments
    filc = n = 0;
    
    // allocate space for max. number of files
    filv = malloc(argc*sizeof(char*));
    if(!filv)
    {
        perror(NULL);
        return EXIT_FAILURE;
    }
    
    // parse arguments
    err = 0;
    for(int i = 1; i < argc && !err; ++i)
    {
        if(argv[i][0] == '-')
        {
            // flags
            for(char* c = &argv[i][1]; *c; ++c)
            {
                // group size
                if(*c == 'n')
                {
                    if(!n && i + 1 < argc)
                        n = atoi(argv[++i]);
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
            filv[filc++] = argv[i];
        }
    }
    
    // make sure at least one input file was given
    if(filc < 1)
        err = 1;
    
    if(err)
    {
        fprintf(stderr, "usage: reg2pts [-n NGROUP] FILE1 [ FILE2 ... ]\n");
        return 1;
    }
    
    // read regions from files
    m = 0;
    for(f = 0; f < filc; ++f)
    {
        fp = fopen(filv[f], "r");
        if(!fp)
        {
            perror(filv[f]);
            return 1;
        }
        
        for(line = 1; fgets(linebuf, sizeof(linebuf), fp); ++line)
        {
            double x = 0, y = 0, dx = 1, dy = 1, rho = 0;
            
            len = strcspn(linebuf, TOK_NL);
            
            if(!linebuf[len] && !feof(fp))
            {
                fprintf(stderr, "%s: line %d: line too long\n", filv[f], line);
                exit(1);
            }
            
            linebuf[len] = 0;
            
            begin = strspn(linebuf, TOK_WS);
            
            if(linebuf[begin] == '#' || linebuf[begin] == '\0')
                continue;
            
            tok = strcspn(linebuf+begin, TOK_DL);
            
            code = 0;
            for(size_t i = 0; i != tok; ++i)
                code += linebuf[begin+i];
            
            switch(code)
            {
                case 625: // global
                case 515: // image
                case 861: // physical
                {
                    continue;
                }
                
                case 261: // fk4
                case 262: // fk5
                case 433: // icrs
                case 824: // galactic
                case 845: // ecliptic
                {
                    // world coordinates
                    goto wcs_error;
                    continue;
                }
                
                case 554: // point
                {
                    if(sscanf(linebuf+begin, "point(%lf,%lf)", &x, &y) != 2)
                        goto syntax_error;
                    break;
                }
                
                case 750: // ellipse
                {
                    double rx, ry, a, c, s;
                    if(sscanf(linebuf+begin, "ellipse(%lf,%lf,%lf,%lf,%lf)", &x, &y, &rx, &ry, &a) != 5)
                        goto syntax_error;
                    c = cos(a*DEGREE);
                    s = sin(a*DEGREE);
                    dx = sqrt(rx*rx*c*c + ry*ry*s*s);
                    dy = sqrt(rx*rx*s*s + ry*ry*c*c);
                    rho = ((rx*rx - ry*ry)*c*s)/(dx*dy);
                    dx *= SCALE;
                    dy *= SCALE;
                    break;
                }
                
                default:
                {
                    // unknown keyword
                    goto keyword_error;
                }
            }
            
            printf("% 18.8f  % 18.8f", x, y);
            if(dx != 1 || dy != dx || rho != 0)
                printf("  % 18.8f", dx);
            if(dy != dx || rho != 0)
                printf("  % 18.8f", dy);
            if(rho != 0)
                printf("  % 18.8f", rho);
            printf("\n");
            
            m += 1;
            
            if(n > 0 && m % n == 0)
                printf("\n");
        }
        
        fclose(fp);
    }
    
    return EXIT_SUCCESS;
    
wcs_error:
    fclose(fp);
    fprintf(stderr, "%s: line %d: world coordinates not supported\n", filv[f], line);
    return EXIT_FAILURE;
    
syntax_error:
    fclose(fp);
    fprintf(stderr, "%s: line %d: syntax error\n", filv[f], line);
    return EXIT_FAILURE;
    
keyword_error:
    fclose(fp);
    fprintf(stderr, "%s: line %d: unknown keyword: %.*s [%u]\n", filv[f], line, (int)tok, linebuf+begin, code);
    return EXIT_FAILURE;
}

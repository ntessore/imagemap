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
    int filc = argc - 1;
    char** filv = argv + 1;
    int file;
    
    FILE* fp;
    int line;
    char linebuf[LINELEN];
    size_t begin, len, tok;
    
    unsigned int code;
    
    if(filc < 1)
    {
        fprintf(stderr, "usage: regpts FILE1 [ FILE2 ... ]\n");
        return 1;
    }
    
    for(file = 0; file < filc; ++file)
    {
        fp = fopen(filv[file], "r");
        if(!fp)
        {
            perror(filv[file]);
            return 1;
        }
        
        for(line = 1; fgets(linebuf, sizeof(linebuf), fp); ++line)
        {
            double x = 0, y = 0, dx = 1, dy = 1, rho = 0;
            
            len = strcspn(linebuf, TOK_NL);
            
            if(!linebuf[len] && !feof(fp))
            {
                fprintf(stderr, "%s: line %d: line too long\n", filv[file], line);
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
                case 861: // physical
                {
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
                    fprintf(stderr, "    # %.*s [%u]\n", (int)tok, linebuf+begin, code);
                    continue;
                }
            }
            
            printf("; %g, %g", x, y);
            if(dx != 1 || dy != 1 || rho != 0)
                printf(", %g", dx);
            if(dy != dx || rho != 0)
                printf(", %g", dy);
            if(rho != 0)
                printf(", %g", rho);
        }
        
        printf("\n");
        
        fclose(fp);
    }
    
    return 0;
    
syntax_error:
    fclose(fp);
    fprintf(stderr, "%s: line %d: syntax error\n", filv[file], line);
    return 1;
}

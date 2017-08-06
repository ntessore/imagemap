#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mpis.h"

// epsilon for stabilising the Cholesky decomposition
#ifndef CHOL_EPS
#define CHOL_EPS 1E-100
#endif

// generate array of random normal variates using the polar method
static void randnv(int c, double v[])
{
    double x[2];
    double s;
    
    for(int i = 0; i < c; ++i)
    {
        if(i%2 == 0)
        {
            do
            {
                x[0] = 2*(rand()/(1. + RAND_MAX)) - 1;
                x[1] = 2*(rand()/(1. + RAND_MAX)) - 1;
                s = x[0]*x[0] + x[1]*x[1];
            }
            while(s >= 1);
        }
        v[i] = x[i%2]*sqrt(-2*log(s)/s);
    }
}

// Cholesky decomposition
static void chol_eps(int n, const double A[], double eps, double L[])
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < i + 1; ++j)
        {
            double s = 0;
            for(int k = 0; k < j; ++k)
                s += L[i*n+k]*L[j*n+k];
            L[i*n+j] = j < i ? (A[i*n+j] - s)/L[j*n+j]
                                     : sqrt(A[i*n+i] + eps - s);
        }
        for(int j = i + 1; j < n; ++j)
            L[i*n+j] = 0;
    }
}

int mpis(mpis_func fn, int m, int n, const double mean[], const double covar[],
         void* private, int num_samples, double samples[])
{
    // sample probability
    double q;
    
    // total probability of samples (= evidence)
    double w;
    
    // status of user function
    int status;
    
    // log-likelihood of sample
    double loglike;
    
    // largest sample probability
    double qmax = -HUGE_VAL;
    
    // allocate arrays
    double* s = malloc(n*n*sizeof(double));
    double* v = malloc(n*sizeof(double));
    double* d = malloc(m*sizeof(double));
    if(!s || !v || !d)
        return MPIS_ERR_MEMORY;
    
    // compute square root of covariance matrix
    chol_eps(n, covar, CHOL_EPS, s);
    
    // draw samples with given mean and variance
    for(int i = 0; i < num_samples; ++i)
    {
        // draw uncorrelated random normal variates
        randnv(n, v);
        
        // compute log-probability of draw, ignore fixed parameters
        q = 0;
        for(int j = 0; j < n; ++j)
            q += -0.5*v[j]*v[j]*(covar[j*n+j] > 0);
        
        // transform to sample using mean and sqrt(covar)
        for(int j = 0; j < n; ++j)
        {
            double x = mean[j];
            for(int k = 0; k < n; ++k)
                x += s[j*n+k]*v[k];
            samples[i*(n+2)+j] = x;
        }
        
        // compute deviates
        status = fn(m, n, &samples[i*(n+2)], d, NULL, private);
        
        // check for user function error
        if(status < 0)
            return status;
        
        // compute log-likelihood of parameters
        loglike = 0;
        for(int j = 0; j < m; ++j)
            loglike += -0.5*d[j]*d[j];
        
        // compute relative sample log-probability
        q = loglike - q;
        
        // store sample likelihood and probability
        samples[i*(n+2)+n+0] = loglike;
        samples[i*(n+2)+n+1] = q;
        
        // store largest sample log-probability
        if(q > qmax)
            qmax = q;
    }
    
    // free arrays
    free(s);
    free(v);
    free(d);
    
    // sum relative sample log-probabilities
    w = 0;
    for(int i = 0; i < num_samples; ++i)
        w += exp(samples[i*(n+2)+n+1] - qmax);
    w = qmax + log(w);
    
    // compute sample probabilities
    for(int i = 0; i < num_samples; ++i)
        samples[i*(n+2)+n+1] = exp(samples[i*(n+2)+n+1] - w);
    
    // return positive value to indicate success
    return 1;
}

double mpis_neff(int n, int num_samples, const double samples[])
{
    // sum of weights and squared weights
    double W, W2;
    
    // sum weights and squared weights
    W = W2 = 0;
    for(int i = 0; i < num_samples; ++i)
    {
        const double w = samples[i*(n+2)+n+1];
        W += w;
        W2 += w*w;
    }
    
    // effective sample size
    return (W*W)/W2;
}

void mpis_stat(int n, int num_samples, const double samples[],
               double mean[], double sigma[])
{
    // cumulative weight
    double cw = 0;
    
    // zero mean and sigma arrays
    memset(mean, 0, n*sizeof(double));
    if(sigma)
        memset(sigma, 0, n*sizeof(double));
    
    // compute mean and variance in one pass
    for(int i = 0; i < num_samples; ++i)
    {
        // weight of sample (= probability)
        const double w = samples[i*(n+2)+n+1];
        
        // skip zero-weight samples
        if(!w)
            continue;
        
        // add to cumulative weight
        cw += w;
        
        // add parameters of sample to mean and variance
        for(int j = 0; j < n; ++j)
        {
            // parameter value
            const double x = samples[i*(n+2)+j];
            
            // difference with current mean
            const double delta = x - mean[j];
            
            // update mean
            mean[j] += (w/cw)*delta;
            
            // update variance (if given)
            if(sigma)
                sigma[j] += w*delta*(x - mean[j]);
        }
    }
    
    // compute sigma (if given)
    if(sigma)
    {
        // effective sample size
        const double neff = mpis_neff(n, num_samples, samples);
        
        // Bessel's correction for effective sample size
        const double b = neff/(neff - 1);
        
        // standard deviation from unbiased variance
        for(int j = 0; j < n; ++j)
            sigma[j] = sqrt(b*sigma[j]/cw);
    }
}

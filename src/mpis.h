#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// shadow MPFIT status definitions
#define MPIS_SUCCESS (1)
#define MPIS_ERR_MEMORY (-20)

// fitting function compatible with MPFIT
typedef int (*mpis_func)(int m, int n, double *x, double *fvec,
                         double **dvec, void *private_data);

// importance sampling from normal distribution with given mean and covariance
int mpis(mpis_func fn, int m, int n, const double mean[], const double covar[],
         void* private, int num_samples, double samples[]);

// effective sample size of weighted sample
double mpis_neff(int n, int num_samples, const double samples[]);

// compute mean and sigma from samples
void mpis_stat(int n, int num_samples, const double samples[],
               double mean[], double sigma[]);

#ifdef __cplusplus
} /* extern "C" */
#endif

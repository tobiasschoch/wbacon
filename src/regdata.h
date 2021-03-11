#ifndef _REGDATA_H
#define _REGDATA_H

// data structure
typedef struct regdata_struct {
    int n;
    int p;
    double *w;          // sampling weight
    double *w_sqrt;     // sqrt of weight
    double *x;          // design matrix (raw and weighted)
    double *wx;
    double *y;          // response vector (raw and weighted)
    double *wy;
} regdata;

// structure of estimates
typedef struct estimate_struct {
    double sigma;       // regression scale
    double *weight;     // weight
    double *resid;      // residuals
    double *beta;       // regression coefficient
    double *dist;       // distance
    double *L;          // Cholesky factor
    double *xty;        // X^Ty
} estimate;
#endif

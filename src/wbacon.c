/* Implementation of the weighted BACON algorithm for multivariate outlier
   detection of Billor et al. (2000), with the extension to allow for weighting
   of Béguin and Hulliger (2008)

   Copyright (C) 2020-2021 Tobias Schoch (e-mail: tobias.schoch@gmail.com)

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, a copy is available at
   https://www.gnu.org/licenses/

   Billor N, Hadi AS, Vellemann PF (2000). BACON: Blocked Adaptative
      Computationally Efficient Outlier Nominators. Computational Statistics
      and Data Analysis 34, pp. 279-298.
   Béguin C, Hulliger B (2008). The BACON-EEM Algorithm for Multivariate
      Outlier Detection in Incomplete Survey Data. Survey Methodology 34,
      pp. 91-103.
*/

#include "wbacon.h"
#define _POWER2(_x) ((_x) * (_x))

// data structure
typedef struct wbdata_struct {
    int n;
    int p;
    double *x;
    double *w;
    double *w_sqrt;
    double *dist;
} wbdata;

// structure of working arrays
typedef struct workarray_struct {
    int *iarray;
    double *work_n;
    double *work_np;
    double *work_pp;
    double *work_2n;
} workarray;

// declarations of local function
static wbacon_error_type initial_subset(wbdata*, workarray*, double* restrict,
    double* restrict, double* restrict, int* restrict, int* restrict, int*,
    int*);
static wbacon_error_type initial_location(wbdata*, workarray*,
    double* restrict, double* restrict, double* restrict, int*);
static inline wbacon_error_type mahalanobis(wbdata*, workarray*,
    double* restrict, double* restrict, double* restrict);
static wbacon_error_type check_matrix_fullrank(double* restrict, int);
static inline void mean_scatter_w(wbdata*, double* restrict, double* restrict,
    double* restrict, double* restrict, double* restrict);
static inline void scatter_w(wbdata*, double* restrict, double* restrict,
    double* restrict, double* restrict);
static inline void euclidean_norm2(wbdata*, double* restrict, double* restrict);
static void verbose_message(int, int, int, double);
static inline double cutoffval(int, int, int) __attribute__((always_inline));

/******************************************************************************\
|* Quantile of chi-square distr. (approximation of Severo and Zelen, 1960)    *|
|*  p    probability                                                          *|
|*  df   degrees of freedom                                                   *|
|*                                                                            *|
|* Severo, N.C. and Zelen M. (1960). Normal Approximation to the Chi-Square   *|
|* and Non-Central F Probability Functions, Biometrika 47, pp. 411-416        *|
\******************************************************************************/
#if 0
static inline double qchisq2(double p, double df)
{
    const double c1 = 0.0740740740;
    const double c2 = 2.8284271247;
    const double c3 = 0.4714045207;
    double xp, xp2, hv, sqrt_df, d_df, res;
    xp = qnorm(p, 0.0, 1.0, 0, 0);
    xp2 = xp * xp;
    d_df = (double)df;
    sqrt_df = sqrt(d_df);
    hv = - c1 / d_df * (c2 * (xp2 - 1.0) / (3.0 * sqrt_df)
        - (xp * xp2 - 3.0 * xp) / 4.0);
    res = 1.0 - 2.0 / (9.0 * d_df) + (xp - hv) * c3 / sqrt_df;
    return d_df * res * res * res;
}
#endif

/******************************************************************************\
|* weighted BACON                                                             *|
|*  x        data, array[n, p]                                                *|
|*  w        weights, array[n]                                                *|
|*  center   on return: array[p]                                              *|
|*  scatter  on return: array[p, p]                                           *|
|*  dist     on return: array[n]                                              *|
|*  n, p     dimensions                                                       *|
|*  alpha    prob.                                                            *|
|*  subset   array[n]                                                         *|
|*  cutoff   on return: chi-square cutoff threshold                           *|
|*  maxiter  on entry: maximal no. of iterations, on return: effective no.    *|
|*  verbose  0: quiet; 1: verbose                                             *|
|*  version2 1: 'Version 2' init. of Billor et al. (2000); 0: 'Version 1'     *|
|*  collect  on entry: parameter to specify the size of the intial subset     *|
|*  success  on return: 1: successful; 0: failure                             *|
|*  threads  set the max number of threads for OpenMP                         *|
\******************************************************************************/
void wbacon(double *x, double *w, double *center, double *scatter, double *dist,
    int *n, int *p, double *alpha, int *subset, double *cutoff, int *maxiter,
    int *verbose, int *version2, int *collect, int *success, int *threads)
{
    int subsetsize, default_no_threads;
    wbacon_error_type err;
    int* restrict subset0 = (int*) Calloc(*n, int);
    double* select_weight = (double*) Calloc(*n, double);

    *success = 1;

    // square root of the weights
    double* w_sqrt = (double*) Calloc(*n, double);
    for (int i = 0; i < *n; i++)
        w_sqrt[i] = sqrt(w[i]);

    // initialize and populate the struct 'wbdata'
    wbdata data;
    wbdata *dat = &data;

    dat->n = *n;
    dat->p = *p;
    dat->x = x;
    dat->w = w;
    dat->w_sqrt = w_sqrt;
    dat->dist = dist;

    // initialize and populate the struct 'workarray'
    workarray warray;
    workarray *work = &warray;

    int *iarray = (int*) Calloc(*n, int);
    double *work_n = (double*) Calloc(*n, double);
    double *work_np = (double*) Calloc(*n * *p, double);
    double *work_pp = (double*) Calloc(*p * *p, double);
    double *work_2n = (double*) Calloc(2 * *n, double);
    work->iarray = iarray;
    work->work_n = work_n;
    work->work_np = work_np;
    work->work_pp = work_pp;
    work->work_2n = work_2n;

    // store current definition of max number of threads
    default_no_threads = omp_get_max_threads();
    // set preferred number of threads
    if (*threads <= default_no_threads) {
        omp_set_num_threads(*threads);
    } else {
        PRINT_OUT("The requested no. of threads is larger than the default.\n");
        PRINT_OUT("Thus, the default is kept at %d\n", default_no_threads);
    }

    // STEP 0
    // initial location
    err = initial_location(dat, work, select_weight, center, scatter, version2);
    if (err != WBACON_ERROR_OK) {
        *success = 0;
        PRINT_OUT("Error: covariance %s\n", wbacon_error(err));
        goto clean_up;
    }

    // initial subset
    err = initial_subset(dat, work, select_weight, center, scatter, subset,
        &subsetsize, verbose, collect);
    if (err != WBACON_ERROR_OK) {
        *success = 0;
        PRINT_OUT("Error: %s (initial subset)\n", wbacon_error(err));
        goto clean_up;
    }

    // STEP 1: update iteratively
    double chi2 = qchisq(*alpha / (double)*n , (double)(*p), 0, 0);
    int iter = 1, is_different;
    for (;;) {
        if (*verbose)
            verbose_message(subsetsize, *n, iter, *cutoff);

        // location, scatter and the Mahalanobis distances
        err = mahalanobis(dat, work, select_weight, center, scatter);
        if (err != WBACON_ERROR_OK) {
            *success = 0;
            PRINT_OUT("Error: covariance %s (iterative updating)\n",
                wbacon_error(err));
            goto clean_up;
        }

        // check whether the subsets differ (XOR current with previous subset)
        is_different = 0;
        for (int i = 0; i < *n; i++) {
            if (subset0[i] ^ subset[i]) {
                is_different = 1;
                break;
            }
        }
        if (is_different == 0) {
            *maxiter = iter;
            break;
        }

        // chi-square cutoff value (quantile)
        *cutoff = chi2 * cutoffval(*n, subsetsize, *p);

        // generate new subset (based on updated Mahalanobis dist.)
        Memcpy(subset0, subset, *n);

        subsetsize = 0;
        for (int i = 0; i < *n; i++) {
            if (dist[i] < *cutoff) {        // obs. is in subset
                subset[i] = 1;
                select_weight[i] = 1.0;
                subsetsize += 1;
        } else {                            // obs. is not in subset
                subset[i] = 0;
                select_weight[i] = 0.0;
            }
        }

        iter++;
        if (iter > *maxiter) {
            *success = 0;
            break;
        }
    }

    for (int i = 0; i < *n; i++)            // Mahalanobis distances
        dist[i] = sqrt(dist[i]);

clean_up:
    Free(subset0); Free(work_np); Free(work_pp);
    Free(work_2n); Free(work_n); Free(iarray); Free(w_sqrt);
    Free(select_weight);

    // set the number of threads to the default value
    if (*threads != default_no_threads)
        omp_set_num_threads(default_no_threads);
}

/******************************************************************************\
|* initial location: either V1 or V2 of Billor et al. (2000)                  *|
|*  dat           data, typedef struct wbdata                                 *|
|*  work          work arrays, typedef struct workarray                       *|
|*  select_weight weight = 1.0 if obs. in subset, otherwise 0.0, array[n]     *|
|*  center        array[p]                                                    *|
|*  scatter       array[p, p]                                                 *|
|*  version2      toogle: 1 = version V2; 0 = Version V1                      *|
\******************************************************************************/
static wbacon_error_type initial_location(wbdata *dat, workarray *work,
    double* restrict select_weight, double* restrict center,
    double* restrict scatter, int* version2)
{
    int n = dat->n, p = dat->p;
    double* restrict x = dat->x;
    wbacon_error_type err = WBACON_ERROR_OK;

    if (*version2) {
        // center: coordinate-wise weighted median
        double d_half = 0.5;
        for (int j = 0; j < p; j++)
            wquantile_noalloc(x + n * j, dat->w, work->work_2n, &n,
                &d_half, &center[j]);

        // distance: Euclidean norm
        euclidean_norm2(dat, work->work_np, center);
    } else {
        // Mahalanobis distances
        for (int i = 0; i < n; i++)
            select_weight[i] = 1.0;

        err = mahalanobis(dat, work, select_weight, center, scatter);
    }
    return err;
}

/******************************************************************************\
|* cutoff for the chi-squared quantile                                        *|
|*  k       size of subset                                                    *|
|*  n, p    dimensions                                                        *|
\******************************************************************************/
static inline double cutoffval(int n, int k, int p)
{
    double dn = (double)n, dp = (double)p, dk = (double)k;
    double h = (dn + dp + 1.0) / 2.0;
    double chr = fmax(0.0, (h - dk) / (h + dk));
    double cnp = 1.0 + (dp + 1.0) / (dn - dp) + 2.0 / (dn - 1.0 - 3.0 * dp);
    return _POWER2(cnp + chr);
}

/******************************************************************************\
|* initial_subset                                                              *|
|*  dat           data, typedef struct wbdata                                 *|
|*  work          work arrays, typedef struct workarray                       *|
|*  select_weight weight = 1.0 if obs. in subset, otherwise 0.0, array[n]     *|
|*  center        array[p]                                                    *|
|*  scatter       array[p, p]                                                 *|
|*  subset        on return: initial subset                                   *|
|*  subsetsize    on return: size initial subset                              *|
|*  verbose       0: quiet; 1: verbose                                        *|
|*  collect       parameter to specify the size of the intial subset          *|
\******************************************************************************/
static wbacon_error_type initial_subset(wbdata *dat, workarray *work,
    double* restrict select_weight, double* restrict center,
    double* restrict scatter, int* restrict subset,
    int* restrict subsetsize, int *verbose, int *collect)
{
    int n = dat->n, p = dat->p;
    int* restrict iarray = work->iarray;
    double* restrict work_np = work->work_np;
    wbacon_error_type status = WBACON_ERROR_OK;

    // determine subset size
    int m = (int)fmin((double)(*collect) * (double)p, (double)n * 0.5);

    // (partially) sort the Mahalanobis distances
    psort_array(dat->dist, iarray, n, m);

    // set weights of observations (m+1):n to zero
    for (int i = 0; i < n; i++)
        select_weight[i] = 0.0;

    // set subset = 1 for the first 1:m elements
    int at;
    for (int i = 0; i < m; i++) {
        at = iarray[i];
        select_weight[at] = 1.0;
        subset[at] = 1;
    }

    // check if scatter matrix has full rank; if not, enlarge the subset
    int next_obs;
    while (m < n) {
        scatter_w(dat, work_np, select_weight, center, scatter);
        status = check_matrix_fullrank(scatter, p);
        if (status == WBACON_ERROR_OK)
            break;

        if (*verbose)
        PRINT_OUT("Initial subset: scatter is rank deficient\n");

        m++;
        next_obs = iarray[m - 1];
        select_weight[next_obs] = 1.0;
        subset[next_obs] = 1;
    }

    *subsetsize = m;
    return status;
}

/******************************************************************************\
|* check if a symmetric pos. def. matrix has full rank p (by Cholesky decomp.)*|
|*  x       array[p, p]                                                       *|
|*  p       dimension                                                         *|
\******************************************************************************/
static wbacon_error_type check_matrix_fullrank(double* restrict x, int p)
{
    int rank = 0;

    // check whether some variances are virtually zero
    for (int i = 0; i < p; i++)
        rank += x[i * (p + 1)] > _RANK_TOLERANCE ? 1 : 0;
    if (rank != p)
        return WBACON_ERROR_NOT_POSITIVE_DEFINITE;

    // Cholesky decomposition
    int info;
    F77_CALL(dpotrf)("L", &p, x, &p, &info);
    if (info != 0)
        return WBACON_ERROR_NOT_POSITIVE_DEFINITE;

    // check whether the diagonal elements of the Cholesky factor are > tol.
    rank = 0;
    for (int i = 0; i < p; i++)
        rank += x[i * (p + 1)] > _RANK_TOLERANCE ? 1 : 0;

    if (rank == p)
        return WBACON_ERROR_OK;
    else
        return WBACON_ERROR_RANK_DEFICIENT;
}

/******************************************************************************\
|* prints out messagges in verbose mode                                       *|
|*  subsetsize   size of the current subset                                   *|
|*  n            number of observations                                       *|
|*  inter        current iteration                                            *|
|*  cutoff       chi-squared quantile/ cutoff value                           *|
\******************************************************************************/
static void verbose_message(int subsetsize, int n, int iter, double cutoff)
{
    double percentage = 100.0 * (double)subsetsize / (double)n;
    if (iter > 1)
        PRINT_OUT("Subset %d: m = %d (%.1f%%); cutoff: %.2f\n", iter,
            subsetsize, percentage, cutoff);
    else
        PRINT_OUT("Subset %d: m = %d (%.1f%%)\n", iter, subsetsize, percentage);
}

/******************************************************************************\
|* weighted covariance/ scatter matrix                                        *|
|*  dat           data, typedef struct wbdata                                 *|
|*  work_np       array[n, p]                                                 *|
|*  select_weight weight = 1.0 if obs. in subset, otherwise 0.0, array[n]     *|
|*  center        array[p]                                                    *|
|*  scatter       on return: array[p, p]                                      *|
\******************************************************************************/
static inline void scatter_w(wbdata *dat, double* restrict work_np,
    double* restrict select_weight, double* restrict center,
    double* restrict scatter)
{
    int n = dat->n, p = dat->p;
    double* restrict x = dat->x;
    double* restrict w = dat->w;
    double* restrict w_sqrt = dat->w_sqrt;

    double sum_w = 0.0;
    for (int i = 0; i < n; i++)
        sum_w += w[i] * select_weight[i];

    // centered data
    #pragma omp parallel for if(n > OMP_MIN_SIZE)
    for (int j = 0; j < p; j++) {
        #pragma omp simd
        for (int i = 0; i < n; i++) {
            work_np[n * j + i] = x[n * j + i] - center[j];
            work_np[n * j + i] *= w_sqrt[i] * select_weight[i];
        }
    }

    // lower triangle of the scatter matrix
    const double d_zero = 0.0;
    double denom = 1.0 / (sum_w - 1.0);
    F77_CALL(dsyrk)("L", "T", &p, &n, &denom, work_np, &n, &d_zero, scatter,
        &p);
}

/******************************************************************************\
|* weighted coordinate-wise mean and covariance/ scatter matrix               *|
|*  dat           data, typedef struct wbdata                                 *|
|*  select_weight weight = 1.0 if obs. in subset, otherwise 0.0, array[n]     *|
|*  work_n        array[n]                                                    *|
|*  work_np       array[n, p]                                                 *|
|*  center        on return: array[p]                                         *|
|*  scatter       on return: array[p, p]                                      *|
\******************************************************************************/
static inline void mean_scatter_w(wbdata *dat, double* restrict select_weight,
    double* restrict work_n, double* restrict work_np, double* restrict center,
    double* restrict scatter)
{
    int n = dat->n, p = dat->p;
    double denom;
    double* restrict x = dat->x;
    double* restrict w = dat->w;
    double* restrict w_sqrt = dat->w_sqrt;

    // sum(w[in subset]) and let work_n[i] = w[i] if in subset, otherwise 0
    double sum_w = 0.0, tmp;
    for (int i = 0; i < n; i++) {
        tmp = select_weight[i] * w[i];
        work_n[i] = tmp;
        sum_w += tmp;
    }

    // coordinate-wise mean and centered data
    denom = 1.0 / sum_w;

    #pragma omp parallel for if(n > OMP_MIN_SIZE)
    for (int j = 0; j < p; j++) {
        center[j] = 0.0;
        #pragma omp simd
        for (int i = 0; i < n; i++)
            center[j] += x[n * j + i] * work_n[i];

        center[j] *= denom;

        // center the data and pre-multiply by sqrt(w[i])
        #pragma omp simd
        for (int i = 0; i < n; i++) {
            work_np[n * j + i] = x[n * j + i] - center[j];
            work_np[n * j + i] *= w_sqrt[i] * select_weight[i];
        }
    }

    // lower triangle of the scatter matrix
    const double d_zero = 0.0;
    denom = 1.0 / (sum_w - 1.0);
    F77_CALL(dsyrk)("L", "T", &p, &n, &denom, work_np, &n, &d_zero, scatter,
        &p);
}

/******************************************************************************\
|* squared Euclidean norm of an array[n, p]                                   *|
|*  dat     data, typedef struct wbdata                                       *|
|*  work_np array[n, p]                                                       *|
|*  center  array[p]                                                          *|
|* NOTE: the algorithm follows S. Hammarling's implementaion of LAPACK:dnorm2 *|
\******************************************************************************/
static inline void euclidean_norm2(wbdata *dat, double* restrict work_np,
    double* restrict center)
{
    int n = dat->n, p = dat->p;
    double abs, scale, ssq;
    double* restrict x = dat->x;
    double* restrict dist = dat->dist;

    for (int i = 0; i < n; i++) {
        ssq = 1.0;
        scale = 0.0;
        for (int j = 0; j < p; j++) {
            work_np[n * j + i] = x[n * j + i] - center[j];
            abs = fabs(work_np[n * j + i]);
            if (abs < DBL_EPSILON)
                continue;
            if (scale <= abs) {
                ssq = 1.0 + ssq * _POWER2(scale / abs);
                scale = abs;
            } else {
                ssq += _POWER2(abs / scale);
            }
        }
        dist[i] = _POWER2(scale) * ssq;
    }
}

/******************************************************************************\
|* Mahalanobis distances                                                      *|
|*  dat           data, typedef struct wbdata                                 *|
|*  work          work array, typedef struct work                             *|
|*  select_weight weight = 1.0 if obs. in subset, otherwise 0.0, array[n]     *|
|*  center        array[p]                                                    *|
|*  scatter       array[p, p]                                                 *|
|* NOTE: on return: dat->dist                                                 *|
\******************************************************************************/
static inline wbacon_error_type mahalanobis(wbdata *dat, workarray *work,
    double* restrict select_weight, double* restrict center,
    double* restrict scatter)
{
    int n = dat->n, p = dat->p;
    double* restrict dist = dat->dist;
    double* restrict x = dat->x;
    double* restrict work_np = work->work_np;

    // coordinate-wise mean and scatter matrix
    mean_scatter_w(dat, select_weight, work->work_n, work_np, center, scatter);

    // Cholesky decomposition of scatter matrix
    Memcpy(work->work_pp, scatter, p * p);
    int info;
    F77_CALL(dpotrf)("L", &p, work->work_pp, &p, &info);
    if (info != 0)
        return WBACON_ERROR_RANK_DEFICIENT;

    // center the data
    #pragma omp parallel for if(n > OMP_MIN_SIZE)
    for (int j = 0; j < p; j++) {
        #pragma omp simd
        for (int i = 0; i < n; i++)
            work_np[n * j + i] = x[n* j + i] - center[j];
    }

    // Solve for y in A * y = B by forward substitution (A = Cholesky factor)
    const double d_one = 1.0;
    F77_CALL(dtrsm)("R", "L", "T", "N", &n, &p, &d_one, work->work_pp, &p,
        work_np, &n);

    // squared Mahalanobis distances (row sums)
    for (int i = 0; i < n; i++)
        dist[i] = _POWER2(work_np[i]);

    for (int j = 1; j < p; j++)
        for (int i = 0; i < n; i++)
            dist[i] += _POWER2(work_np[n * j + i]);

    return WBACON_ERROR_OK;
}
#undef _POWER2

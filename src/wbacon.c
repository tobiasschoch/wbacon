/* Implementation of the weighted BACON algorithm for multivariate outlier
   detection of Billor et al. (2000), with the extension to allow for weighting
   of Béguin and Hulliger (2008)

   Copyright (C) 2020 Tobias Schoch (e-mail: tobias.schoch@gmail.com)

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
      Computationally efficient Outlier Nominators. Computational Statistics
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
wbacon_error_type initialsubset(wbdata*, workarray*, double*, double*, int*,
	int*, int*, int*);
wbacon_error_type mahalanobis(wbdata*, double*, double*, double*, double*);
wbacon_error_type check_matrix_fullrank(double*, int, int);
void weightedmean(wbdata*, double*);
void weightedscatter(wbdata*, double*, double*, double*);
void euclidean_norm2(wbdata*, double*, double*);
void verbose_message(int, int, int, double);
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
\******************************************************************************/
void wbacon(double *x, double *w, double *center, double *scatter, double *dist,
	int *n, int *p, double *alpha, int *subset, double *cutoff,
	int *maxiter, int *verbose, int *version2, int *collect, int *success)
{
	*success = 1;
	wbacon_error_type err;

	int *subset0 = (int*) Calloc(*n, int);
	double *original_w = (double*) Calloc(*n, double);
	Memcpy(original_w, w, *n);

	// initialize and populate the struct 'wbdata'
	wbdata data;
	wbdata *dat = &data;

	dat->n = *n;
	dat->p = *p;
	dat->x = x;
	dat->w = w;
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

	// STEP 0: establish initial subset
	if (*version2) {
		// center: coordinate-wise weighted median
		double d_half = 0.5;
		for (int j = 0; j < *p; j++)
			wquantile_noalloc(x + *n * j, w, work_2n, n, &d_half, &center[j]);

		// distance: Euclidean norm
		euclidean_norm2(dat, work_np, center);
	} else {
		// Mahalanobis distances
		err = mahalanobis(dat, work_np, work_pp, center, scatter);
		if (err != WBACON_ERROR_OK) {
			*success = 0;
			PRINT_OUT("Error: covariance %s\n", wbacon_error(err));
			goto clean_up;
		}
	}

	// determine initial subset
	int subsetsize;
	err = initialsubset(dat, work, center, scatter, subset, &subsetsize,
		verbose, collect);

	if (err != WBACON_ERROR_OK) {
		*success = 0;
		PRINT_OUT("Error: %s (initial subset)\n", wbacon_error(err));
		goto clean_up;
	}

	// STEP 1: update iteratively
	double chi2 = sqrt(qchisq(*alpha / (double)*n , (double)(*p), 0, 0));
	int iter = 1, is_different;
	for (;;) {
		if (*verbose)
			verbose_message(subsetsize, *n, iter, *cutoff);

		// location, scatter and the Mahalanobis distances
		err = mahalanobis(dat, work_np, work_pp, center, scatter);
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
		Memcpy(w, original_w, *n);

		subsetsize = 0;
		for (int i = 0; i < *n; i++) {
			if (dist[i] < *cutoff) {
				subset[i] = 1;
				subsetsize += 1;
			} else {
				subset[i] = 0;
				w[i] = 0.0;			// weight = 0 (if obs. is not in subset)
			}
		}

		iter++;
		if (iter > *maxiter) {
			*success = 0;
			break;
		}
	}

clean_up:
	Free(subset0); Free(original_w); Free(work_np); Free(work_pp);
	Free(work_2n); Free(work_n); Free(iarray);
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
	return cnp + chr;
}

/******************************************************************************\
|* initialsubset                                                              *|
|*  dat        data, typedef struct wbdata                                    *|
|*  work       work arrays, typedef struct workarray                          *|
|*  center     array[p]                                                       *|
|*  scatter    array[p, p]                                                    *|
|*  subset     on return: initial subset                                      *|
|*  subsetsize on return: size initial subset                                 *|
|*  verbose    0: quiet; 1: verbose                                           *|
|*  collect    parameter to specify the size of the intial subset             *|
\******************************************************************************/
wbacon_error_type initialsubset(wbdata *dat, workarray *work, double *center,
	double *scatter, int *subset, int *subsetsize, int *verbose, int *collect)
{
	int n = dat->n, p = dat->p;
	wbacon_error_type status = WBACON_ERROR_OK;

	// store 'w' (it will be modifed)
	Memcpy(work->work_n, dat->w, n);

	// sort the Mahalanobis distances
	psort_array(dat->dist, work->iarray, n, n);

	// determine subset size
	int m = (int)fmin((double)(*collect) * (double)p, (double)n * 0.5);

	// set weights of observations (m+1):n to zero
	for (int i = m; i < n; i++)
		dat->w[work->iarray[i]] = 0.0;

	// set subset = 1 for the first 1:m elements
	for (int i = 0; i < m; i++)
		subset[work->iarray[i]] = 1;

	// check if scatter matrix has full rank; if not, enlarge the subset
	while (m < n) {
		weightedscatter(dat, work->work_np, center, scatter);
		status = check_matrix_fullrank(scatter, p, 1);
		if (status == WBACON_ERROR_OK)
			break;

		if (*verbose)
			PRINT_OUT("Initial subset: scatter is rank deficient\n");

		m++;
		int next_obs = work->iarray[m - 1];
		dat->w[next_obs] = work->work_n[next_obs];
		subset[next_obs] = 1;
	}

	*subsetsize = m;
	return status;
}

/******************************************************************************\
|* check if a symmetric pos. def. matrix has full rank p (by Cholesky decomp.)*|
|*  x       array[p, p]                                                       *|
|*  p       dimension                                                         *|
|*  decomp  1: chol. decomposition is computed; 0: x is the chol. decomp.     *|
\******************************************************************************/
wbacon_error_type check_matrix_fullrank(double *x, int p, int decom)
{
	int rank = 0;
	if (decom) {
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
	}

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
void verbose_message(int subsetsize, int n, int iter, double cutoff)
{
	double percentage = 100.0 * (double)subsetsize / (double)n;
	if (iter > 1)
		PRINT_OUT("Subset %d: n = %d (%.1f%%); cutoff: %.2f\n", iter,
			subsetsize, percentage, cutoff);
	else
		PRINT_OUT("Subset %d: n = %d (%.1f%%)\n", iter, subsetsize, percentage);
}

/******************************************************************************\
|* weighted mean (vector valued)                                              *|
|*  dat     data, typedef struct wbdata                                       *|
|*  center  on return: array[p]                                               *|
\******************************************************************************/
void weightedmean(wbdata *dat, double *center)
{
	for (int i = 0; i < dat->p; i++)
		center[i] = 0.0;

	double sum_w = 0.0;
	for (int i = 0; i < dat->n; i++) {
		sum_w += dat->w[i];
		for (int j = 0; j < dat->p; j++)
			center[j] += dat->x[dat->n * j + i] * dat->w[i];
   	}

	for (int j = 0; j < dat->p; j++)
		center[j] /= sum_w;
}

/******************************************************************************\
|* weighted covariance/ scatter matrix                                        *|
|*  dat     data, typedef struct wbdata                                       *|
|*  work    array[n, p]                                                       *|
|*  center  array[p]                                                          *|
|*  scatter on return: array[p, p]                                            *|
\******************************************************************************/
void weightedscatter(wbdata *dat, double *work, double *center, double *scatter)
{
	int n = dat->n, p = dat->p;
	Memcpy(work, dat->x, n * p);

	// center the data and multiply by sqrt(w)
	double sum_w = 0.0;
	for (int i = 0; i < n; i++) {
		sum_w += dat->w[i];
		for (int j = 0; j < p; j++) {
			work[n * j + i] -= center[j];
			work[n * j + i] *= sqrt(dat->w[i]);
		}
	}

	// cross product: scatter matrix = t(work) %*% work;
	const double d_one = 1.0, d_zero = 0.0;
	F77_CALL(dgemm)("T", "N", &p, &p, &n, &d_one, work, &n, work, &n, &d_zero,
		scatter, &p);

	for (int i = 0; i < (p * p); i++)
		scatter[i] /= (sum_w - 1.0);
}

/******************************************************************************\
|* squared Euclidean norm of an array[n, p]                                   *|
|*  dat     data, typedef struct wbdata                                       *|
|*  work    array[n, p]                                                       *|
|*  center  array[p]                                                          *|
\******************************************************************************/
void euclidean_norm2(wbdata *dat, double *work, double *center)
{
	int n = dat->n, p = dat->p;
	Memcpy(work, dat->x, n * p);

	for (int i = 0; i < n; i++) {
		dat->dist[i] = 0.0;
		for (int j = 0; j < p; j++) {
			work[n * j + i] -= center[j];
			dat->dist[i] += _POWER2(work[n * j + i]);
		}
	}
}

/******************************************************************************\
|* Mahalanobis distances                                                      *|
|*  dat      data, typedef struct wbdata                                      *|
|*  work_np  work array[n, p]                                                 *|
|*  work_pp  work array[p, p]                                                 *|
|*  center   array[p]                                                         *|
|*  scatter  array[p, p]                                                      *|
\******************************************************************************/
wbacon_error_type mahalanobis(wbdata *dat, double *work_np, double *work_pp,
	double *center, double *scatter)
{
	int n = dat->n, p = dat->p;

	// center and scatter
	weightedmean(dat, center);
	weightedscatter(dat, work_np, center, scatter);

	Memcpy(work_np, dat->x, n * p);					// copy of 'x'

	for (int i = 0; i < n; i++)
		for (int j = 0; j < p; j++)
			work_np[n * j + i] -= center[j];		// center the data

	// Cholesky decomposition of scatter matrix
	Memcpy(work_pp, scatter, p * p);

	int info;
	F77_CALL(dpotrf)("L", &p, work_pp, &p, &info);
	if (info != 0)
		return WBACON_ERROR_RANK_DEFICIENT;

	// Solve for y in A * y = B by forward substitution (A = Cholesky factor)
	const double d_one = 1.0;
	F77_CALL(dtrsm)("R", "L", "T", "N", &n, &p, &d_one, work_pp, &p, work_np,
		&n);

	// Mahalanobis distances (row sums)
	for (int i = 0; i < n; i++) {
		dat->dist[i] = 0.0;
		for (int j = 0; j < p; j++)
			dat->dist[i] += _POWER2(work_np[n * j + i]);

		dat->dist[i] = sqrt(dat->dist[i]);
	}

	return WBACON_ERROR_OK;
}
#undef _POWER2

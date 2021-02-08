/* Implementation of the weighted BACON algorithm of multivariate outlier
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

// declarations of local function
wbacon_error_type initialsubset(double*, double*, double*, int*, int*, int*, 
	int*, int*);
wbacon_error_type mahalanobis(double*, double*, double*, double*, double*, 
	double*, int*, int*);
wbacon_error_type check_matrix_fullrank(double*, int*, int);
void weightedmean(double*, double*, double*, int*, int*);
void weightedscatter(double*, double*, double*, double*, double*, int*, int*);
void euclidean_norm2(double*, double*, double*, double*, int*, int*);
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
|*  x          array[n, p]                                                    *|
|*  w          on return: elements not in subset have w = 0, array[n]         *|
|*  dist       array[n]                                                       *|
|*  subset     on return: initial subset                                      *|
|*  subsetsize on return: size initial subset                                 *|
|*  n, p       dimensions                                                     *|
|*  verbose    0: quiet; 1: verbose                                           *|
\******************************************************************************/
wbacon_error_type initialsubset(double *x, double *w, double *dist, int *subset, 
	int *subsetsize, int *n, int *p, int *verbose)
{
	int *iarray;
	double *dist_sorted, *x_sorted, *w_sorted, *w_sortedcpy, *center, *scatter,
		*work;
	wbacon_error_type status = WBACON_ERROR_OK;

	center = (double*) Calloc(*p, double);
	scatter = (double*) Calloc(*p * *p, double);
	work = (double*) Calloc(*n * *p, double);

	// STEP 0: sort Mahalanobis dist. (iarray is sorted along with it)
	dist_sorted = (double*) Calloc(*n, double);
	Memcpy(dist_sorted, dist, *n);
	iarray = (int*) Calloc(*n, int);
	psort_array(dist_sorted, iarray, *n, *n);

	// use iarray to sort x (x_sorted) and w (w_sorted)
	x_sorted = (double*) Calloc(*n * *p, double);
	w_sorted = (double*) Calloc(*n, double);
	Memcpy(x_sorted, x, *n * *p);
	Memcpy(w_sorted, w, *n);

	for (int i = 0; i < *n; i++) {
		w_sorted[i] = w[iarray[i]];
		for (int j = 0; j < *p; j++)
			x_sorted[i + *n * j] = x[iarray[i] + *n * j];
	}

	// determine subset size
	int m = (int)fmin(4.0 * (double)*p, (double)*n * 0.5);

	w_sortedcpy = (double*) Calloc(*n, double);
	Memcpy(w_sortedcpy, w_sorted, *n);

	// STEP 1 check if scatter matrix (with obs. 1:m) has full rank
	while (m < *n) {
		// set weights of observations (m+1):n to zero 
		for (int i = m; i < *n; i++)
			w_sortedcpy[i] = 0.0;

		weightedscatter(x_sorted, work, w_sortedcpy, center, scatter, n, p);
		status = check_matrix_fullrank(scatter, p, 1);
		if (status == WBACON_ERROR_OK)		
			break;
		m++;
		Memcpy(w_sortedcpy, w_sorted, *n);
	}

	// STEP 2: generate initial subset and set w = 0 if not in subset
	for (int i = 0; i < m; i++)
		subset[iarray[i]] = 1;
   
	for (int i = m; i < *n; i++)	 
		w[iarray[i]] = 0.0;

	*subsetsize = m;

	Free(iarray); Free(x_sorted); Free(w_sorted); Free(w_sortedcpy);
	Free(center); Free(dist_sorted); Free(scatter); Free(work);

	return status;
}

/******************************************************************************\
|* check if a symmetric pos. def. matrix has full rank p (by Cholesky decomp.)*|
|*  x       array[p, p]                                                       *|
|*  p       dimension                                                         *|
|*  decomp  1: chol. decomposition is computed; 0: x is the chol. decomp.     *|
\******************************************************************************/
wbacon_error_type check_matrix_fullrank(double *x, int *p, int decom)
{
	int rank = 0;
	if (decom) {
		// check whether some variances are virtually zero
		for (int i = 0; i < *p; i++)
			rank += x[i * (*p + 1)] > _RANK_TOLERANCE ? 1 : 0;
		if (rank != *p)
			return WBACON_ERROR_NOT_POSITIVE_DEFINITE;

		// Cholesky decomposition 
		int info;
		F77_CALL(dpotrf)("L", p, x, p, &info);
		if (info != 0) 
			return WBACON_ERROR_NOT_POSITIVE_DEFINITE;
	}

	// check whether the diagonal elements of the Cholesky factor are > tol.
	rank = 0;
	for (int i = 0; i < *p; i++)
		rank += x[i * (*p + 1)] > _RANK_TOLERANCE ? 1 : 0;

	if (rank == *p)
		return WBACON_ERROR_OK;
	else 
		return WBACON_ERROR_RANK_DEFICIENT; 
}

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
|*  success   on return: 1: successful; 0: failure                            *|
\******************************************************************************/
void wbacon(double *x, double *w, double *center, double *scatter, double *dist,
	int *n, int *p, double *alpha, int *subset, double *cutoff,
	int *maxiter, int *verbose, int *version2, int *success)
{
	int *subset0;
	double chi2;
	double *w_cpy, *work_np, *work_pp, *work_2n;
	wbacon_error_type err;
	*success = 1;

	subset0 = (int*) Calloc(*n, int);
	w_cpy = (double*) Calloc(*n, double);
	work_np = (double*) Calloc(*n * *p, double);
	work_pp = (double*) Calloc(*p * *p, double);
	work_2n = (double*) Calloc(2 * *n, double);

	// STEP 0: establish initial subset 
	if (*version2) {
		// center: coordinate-wise weighted median
		double d_half = 0.5;
		for (int j = 0; j < *p; j++)	 
			wquantile_noalloc(x + *n * j, w, work_2n, n, &d_half, &center[j]);

		// distance: Euclidean norm
		euclidean_norm2(x, dist, work_np, center, n, p);

	} else {
		// center: weighted mean
		weightedmean(x, w, center, n, p);

		// scatter and Mahalanobis distances
		weightedscatter(x, work_np, w, center, scatter, n, p);
		err = mahalanobis(x, work_np, work_pp, center, scatter, dist, n, p);
		if (err != WBACON_ERROR_OK) {
			*success = 0;
			PRINT_OUT("Error: covariance %s\n", wbacon_error(err));
			goto clean_up;
		}
	}

	// copy w (w_cpy will be set 0 if obs. is not in the subset)
	Memcpy(w_cpy, w, *n);  

	// determine initial subset
	int subsetsize;
	err = initialsubset(x, w_cpy, dist, subset, &subsetsize, n, p, verbose);
	if (err != WBACON_ERROR_OK) {
		*success = 0;
		PRINT_OUT("Error: %s (initial subset)\n", wbacon_error(err));
		goto clean_up;
	}

	// STEP 1: update iteratively
	chi2 = sqrt(qchisq(*alpha / (double)*n , (double)(*p), 0, 0));	
	int iter = 1, is_different;
	for (;;) {
		if (*verbose) {
			double percentage = 100.0 * (double)subsetsize / (double)*n;
			if (iter > 1)
				PRINT_OUT("Subset %d: n = %d (%.1f%%); cutoff: %.2f\n", iter, 
					subsetsize, percentage, *cutoff);
			else 
				PRINT_OUT("Subset %d: n = %d (%.1f%%)\n", iter, subsetsize, 
					percentage);
		}

		// location, scatter and the Mahalanobis distances
		weightedmean(x, w_cpy, center, n, p);
		weightedscatter(x, work_np, w_cpy, center, scatter, n, p);
		err = mahalanobis(x, work_np, work_pp, center, scatter, dist, n, p);
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
		Memcpy(w_cpy, w, *n);

		subsetsize = 0;
		for (int i = 0; i < *n; i++) {
			if (dist[i] < *cutoff) {
				subset[i] = 1;
				subsetsize += 1;
			} else {
				subset[i] = 0;
				w_cpy[i] = 0.0;	    // weight = 0 (if obs. is not in subset)
			}
		}

		iter++;
		if (iter > *maxiter) {
			*success = 0;
			break;
		}
	}

clean_up:
	Free(subset0); Free(w_cpy); Free(work_np); Free(work_pp); Free(work_2n);
}

/******************************************************************************\
|* weighted mean (vector values)                                              *|
|*  x       array[n, p]                                                       *|
|*  w       array[n], weights                                                 *|
|*  center  on return: array[p]                                               *|
|*  n, p    dimensions                                                        *|
\******************************************************************************/
void weightedmean(double *x, double *w, double *center, int *n, int *p)
{
	for (int i = 0; i < *p; i++)
		center[i] = 0.0;

	double sum_w = 0.0;
	for (int i = 0; i < *n; i++) {
		sum_w += w[i];
		for (int j = 0; j < *p; j++) 
			center[j] += x[*n * j + i] * w[i];
   	}

	for (int j = 0; j < *p; j++)
		center[j] /= sum_w;
}

/******************************************************************************\
|* weighted covariance/ scatter matrix                                        *|
|*  x       array[n, p]                                                       *|
|*  work    array[n, p]                                                       *|
|*  w       array[n], weights                                                 *|
|*  center  array[p]                                                          *|
|*  scatter on return: array[p, p]                                            *|
|*  n, p    dimensions                                                        *|
\******************************************************************************/
void weightedscatter(double *x, double *work, double *w, double *center,
	double *scatter, int *n, int *p)
{
	Memcpy(work, x, *n * *p);

	// center the data and multiply by sqrt(w)
	double sum_w = 0.0;
	for (int i = 0; i < *n; i++) {
		sum_w += w[i];
		for (int j = 0; j < *p; j++) {
			work[*n * j + i] -= center[j];
			work[*n * j + i] *= sqrt(w[i]);
		}
	}

	// cross product: scatter matrix = t(work) %*% work;
	const double d_one = 1.0, d_zero = 0.0;
	F77_CALL(dgemm)("T", "N", p, p, n, &d_one, work, n, work, n, &d_zero,
		scatter, p);

	for (int i = 0; i < (*p * *p); i++)
		scatter[i] /= (sum_w - 1);
}

/******************************************************************************\
|* squared Euclidean norm of an array[n, p]                                   *|
|*  x       array[n, p]                                                       *|
|*  dist    array[n]                                                          *|
|*  work    array[n, p]                                                       *|
|*  center  array[p]                                                          *|
|*  n, p    dimensions                                                        *|
\******************************************************************************/
void euclidean_norm2(double *x, double *dist, double *work, double *center, 
	int *n, int *p)
{
	Memcpy(work, x, *n * *p);

	for (int i = 0; i < *n; i++) {
		dist[i] = 0.0;
		for (int j = 0; j < *p; j++) {
			work[*n * j + i] -= center[j];
			dist[i] += work[*n * j + i] * work[*n * j + i];
		}
	}
}

/******************************************************************************\
|* Mahalanobis distances                                                      *|
|*  x        array[n, p], n obs. on p. variables, whose M. dist. is computed  *|
|*  work_np  work array[n, p]                                                 *|
|*  work_pp  work array[p, p]                                                 *|
|*  center   array[p]                                                         *|
|*  scatter  array[p, p]                                                      *|
|*  dist     on return: array[n]                                              *|
|*  n, p     dimensions                                                       *|
\******************************************************************************/
wbacon_error_type mahalanobis(double *x, double *work_np, double *work_pp, 
	double *center, double *scatter, double *dist, int *n, int *p)
{
	Memcpy(work_np, x, *n * *p);					// copy of 'x' 

	for (int i = 0; i < *n; i++)
		for (int j = 0; j < *p; j++)
			work_np[*n * j + i] -= center[j];		// center the data

	// Cholesky decomposition of scatter matrix
	Memcpy(work_pp, scatter, *p * *p);
	int info;
	F77_CALL(dpotrf)("L", p, work_pp, p, &info);
	if (info != 0) 
		return WBACON_ERROR_RANK_DEFICIENT;

	// Solve for y in A * y = B by forward substitution (A = Cholesky factor)
	const double d_one = 1.0;
	F77_CALL(dtrsm)("R", "L", "T", "N", n, p, &d_one, work_pp, p, work_np, n);

	// Mahalanobis distances (row sums)
	for (int i = 0; i < *n; i++) {   
		dist[i] = 0.0;
		for (int j = 0; j < *p; j++)
			dist[i] += _POWER2(work_np[*n * j + i]);

		dist[i] = sqrt(dist[i]);
	}

	return WBACON_ERROR_OK; 
}
#undef _POWER2 

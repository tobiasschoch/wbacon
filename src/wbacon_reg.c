/* Implementation of the weighted BACON algorithm for robust linear regression
   of Billor et al. (2000)

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
*/

//FIXME
// regression_scale -> into compute_ti

#include "wbacon_reg.h"
#define _POWER2(_x) ((_x) * (_x))

// structure of working arrays
typedef struct workarray_struct {
	int lwork;
	int *iarray;
	double *work_p;
	double *work_n;
	double *work_np;
	double *work_pp;
	double *dgels_work;
} workarray;

// structure of estimates
typedef struct estimate_struct {
	double sigma;	// regression scale
	double *weight;	// weight
	double *resid;	// residuals
	double *beta;	// regression coefficient
	double *dist;	// distance
	double *L;		// Cholesky factor
	double *xty;	// X^Ty
} estimate;

// declarations of local function
wbacon_error_type initial_reg(regdata*, workarray*, estimate*, int*, int*,
	int*);
wbacon_error_type algorithm_4(regdata*, workarray*, estimate*, int*, int*, int*,
	int*, int*);
wbacon_error_type algorithm_5(regdata *dat, workarray *work, estimate *est,
	int *subset0, int *subset1, double *alpha, int *m, int *maxiter,
	int *verbose);
wbacon_error_type compute_ti(regdata*, workarray*, estimate*, int*, int*,
	double*);
wbacon_error_type update_chol_xty(regdata*, workarray*, estimate*, int*, int*,
	int*);
wbacon_error_type chol_downdate(double*, double*, int);
wbacon_error_type hat_matrix(regdata*, workarray*, double*, double*);
double regression_scale(double*, double*, int, int);
void select_subset(double*, int*, int*, int*, int*);
void cholesky_reg(double*, double*, double*, double*, int*, int*);
void chol_update(double*, double*, int);

/******************************************************************************\
|* BACON regression estimator                                                 *|
\******************************************************************************/
//FIXME
void wbacon_reg(double *x, double *y, double *w, double *resid, double *beta,
	int *subset0, double *dist, int *n, int *p, int *m, int *verbose,
	int *success, int *collect, double *alpha, int *maxiter)
{
	wbacon_error_type err;
	*success = 1;

	int *subset1 = (int*) Calloc(*n, int);

	// initialize and populate 'data' which is a regdata struct
	regdata data;
	regdata *dat = &data;
	dat->n = *n; dat->p = *p;
	dat->x = x;
	dat->y = y;
	dat->w = w;
	double *wy = (double*) Calloc(*n, double);
	dat->wy = wy;
	double *wx = (double*) Calloc(*n * *p, double);
	dat->wx = wx;

	// initialize and populate 'est' which is a estimate struct
	estimate the_estimate;
	estimate *est = &the_estimate;
	double *weight = (double*) Calloc(*n, double);
	est->weight = weight;
	est->resid = resid;
	est->beta = beta;
	est->dist = dist;
	double *L = (double*) Calloc(*p * *p, double);
	est->L = L;
	double *xty = (double*) Calloc(*p, double);
	est->xty = xty;

	// initialize and populate 'data' which is a workarray struct
	workarray warray;
	workarray *work = &warray;
	double *work_p = (double*) Calloc(*p, double);
	work->work_p = work_p;
	double *work_pp = (double*) Calloc(*p * *p, double);
	work->work_pp = work_pp;
	double *work_np = (double*) Calloc(*n * *p, double);
	work->work_np = work_np;
	double *work_n = (double*) Calloc(*n, double);
	work->work_n = work_n;
	int *iarray = (int*) Calloc(*n, int);
	work->iarray = iarray;
	// determine size of work array for LAPACK:degels
	work->lwork = fitwls(dat, w, work_np, -1, est->beta, est->resid);
	double *dgels_work = (double*) Calloc(work->lwork, double);
	work->dgels_work = dgels_work;

	// STEP 0 (initialization)
	err = initial_reg(dat, work, est, subset0, m, verbose);
	if (err != WBACON_ERROR_OK) {
		*success = 0;
		PRINT_OUT("Error: design %s (step 0)\n", wbacon_error(err));
		goto clean_up;
	}

	// compute t[i]'s
	err = compute_ti(dat, work, est, subset0, m, est->dist);
	if (err != WBACON_ERROR_OK) {
		*success = 0;
		PRINT_OUT("Error: %s (step 0)\n", wbacon_error(err));
		goto clean_up;
	}

	// select the p+1 obs. with the smallest ti's (initial basic subset)
	*m = *p + 1;
	select_subset(est->dist, iarray, subset1, m, n);

	// STEP 1 (Algorithm 4)
	err = algorithm_4(dat, work, est, subset0, subset1, m, verbose, collect);
	if (err != WBACON_ERROR_OK) {
		PRINT_OUT("Error: %s (Cholesky update, step 1)\n", wbacon_error(err));
		*success = 0;
		goto clean_up;
	}

	// STEP 2 (Algorithm 5)
	err = algorithm_5(dat, work, est, subset0, subset1, alpha, m, maxiter,
		verbose);
	if (err != WBACON_ERROR_OK) {
		PRINT_OUT("Error: %s (step 2)\n", wbacon_error(err));
		*success = 0;
	}

	// copy the QR factorization to x (as returned by fitwls -> dgels -> dgeqrf)
	Memcpy(x, wx, *n * *p);

clean_up:
	Free(work_pp); Free(work_p); Free(work_np); Free(work_n); Free(weight);
	Free(iarray);
	Free(wx); Free(wy);
	Free(L); Free(subset1); Free(xty); Free(dgels_work);
}

/******************************************************************************\
|* Initial basic subset, adapted for weighting                                *|
|*  dat      typedef struct regdata                                           *|
|*  work     typedef struct workarray                                         *|
|*  est      typedef struct estimate                                          *|
|*  subset   subset, 1: obs. is in the subset, 0 otherwise, array[n]          *|
|*  m        number of obs. in subset                                         *|
|*  verbose  toggle: 1: additional information is printed to the console;     *|
|*           0: quiet                                                         *|
|* On return, dat->wx is overwritten with the R matrix of the QR factorization*|
\******************************************************************************/
wbacon_error_type initial_reg(regdata *dat, workarray *work, estimate *est,
	int *subset, int *m, int *verbose)
{
	int info, n = dat->n, p = dat->p;
	wbacon_error_type status = WBACON_ERROR_OK;

	// 3) set weight = w if obs. is in subset (otherwise weight = 0)
	Memcpy(est->weight, dat->w, n);
	for (int i = 0; i < n; i++)
		if (subset[i] == 0)
			est->weight[i] = 0.0;

	// 4) compute regression estimate (on return, dat->wx is overwritten by the
	// R matrix of the QR factorization (R will be used by the caller of
	// initial_reg)
	info = fitwls(dat, est->weight, work->dgels_work, work->lwork, est->beta,
		est->resid);

	// if the design matrix dat->x on the subset is rank deficient, we enlarge
	// the subset
	if (info) {
		status = WBACON_ERROR_RANK_DEFICIENT;
		// sort the dist[i]'s in ascending order
		psort_array(est->dist, work->iarray, n, n);

		// add obs. to the subset until x has full rank
		while (*m < n) {
			(*m)++;
			// select obs. with smallest dist (among those obs. not in the
			// subset); and set its weight to the original weight (i.e., not 0)
			int at = work->iarray[*m - 1];
			subset[at] = 1;
			est->weight[at] = dat->w[at];

			// re-do regression and check rank
			info = fitwls(dat, est->weight, work->dgels_work, work->lwork,
				est->beta, est->resid);
			if (info == 0) {
				status = WBACON_ERROR_OK;
				break;
			}
		}
	}
	if (*verbose)
		PRINT_OUT("Step 0: initial subset, m = %d\n", *m);

	// extract R matrix (as a lower triangular matrix: L) from dat->wx
	for (int i = 0; i < p; i++)
		for (int j = i; j < p; j++)
			est->L[j + i * p] = dat->wx[i + j * n];

	// compute xty (weighted)
	for (int i = 0; i < p; i++)
		for (int j = 0; j < n; j++)
			if (subset[j])
				est->xty[i] += dat->w[j] * dat->x[j + i * n] * dat->y[j];

	// compute the regression scale (weighted)
	est->sigma = regression_scale(est->resid, est->weight, n, p);

	return status;
}

/******************************************************************************\
|* Algorithm 4 of Billor et al. (2000), adapted for weighting                 *|
|*  dat      typedef struct regdata                                           *|
|*  work     typedef struct workarray                                         *|
|*  est      typedef struct estimate                                          *|
|*  subset0  subset, 1: obs. is in the subset, 0 otherwise, array[n]          *|
|*  subset1  subset, 1: obs. is in the subset, 0 otherwise, array[n]          *|
|*  m        number of obs. in subset1                                        *|
|*  verbose  toggle: 1: additional information is printed to the console;     *|
|*           0: quiet                                                         *|
|*  collect  defines size of initial subset: m = collect * p                  *|
\******************************************************************************/
wbacon_error_type algorithm_4(regdata *dat, workarray *work, estimate *est,
	int *subset0, int *subset1, int *m, int *verbose, int *collect)
{
	int n = dat->n, p = dat->p;
	wbacon_error_type err;

	if (*verbose)
		PRINT_OUT("Step 1 (Algorithm 4):\n");

	// STEP 1 (Algorithm 4)
	while (*m <= p * *collect) {
		if (*verbose)
			PRINT_OUT("  m = %d", *m);

		// update cholesky factor and xty matrix (subset0 => subset1)
		err = update_chol_xty(dat, work, est, subset0, subset1, verbose);

		// check whether L is well defined; if not, keep adding obs. to the
		// subset until it has full rank
		if (err != WBACON_ERROR_OK) {
			for (;;) {
				(*m)++;
				subset1[work->iarray[*m - 1]] = 1;

				if (*verbose)
					PRINT_OUT("  m = %d", *m);

				// re-do the updating
				err = update_chol_xty(dat, work, est, subset0, subset1, verbose);
				if (err == WBACON_ERROR_OK)
					break;
				if (*m == p * 4)
					return err;
			}
		}

		// prepare next while iteration; then update subset1
		Memcpy(subset0, subset1, n);

		// regression estimate: beta (updated Cholesky factor)
		cholesky_reg(est->L, dat->x, est->xty, est->beta, &n, &p);

		// compute residuals
		const int int_1 = 1;
		const double double_minus1 = -1.0, double_1 = 1.0;
		Memcpy(est->resid, dat->y, n);
		F77_CALL(dgemv)("N", &n, &p, &double_minus1, dat->x, &n, est->beta,
			&int_1, &double_1, est->resid, &int_1);

		// compute regression scale (sigma)
		est->sigma = regression_scale(est->resid, est->weight, n, p);

		// compute t[i]'s
		err = compute_ti(dat, work, est, subset1, m, est->dist);
		if (err != WBACON_ERROR_OK)
			return err;

		// select the m obs. with the smallest t[i]'s
		(*m)++;
		select_subset(est->dist, work->iarray, subset1, m, &n);
	}
	return WBACON_ERROR_OK;
}

/******************************************************************************\
|* Algorithm 5 of Billor et al. (2000), adapted for weighting                 *|
\******************************************************************************/
wbacon_error_type algorithm_5(regdata *dat, workarray *work, estimate *est,
	int *subset0, int *subset1, double *alpha, int *m, int *maxiter,
	int *verbose)
{
	int n = dat->n, p = dat->p, iter = 1, info, i;
	double cutoff;
	wbacon_error_type err;

	if (*verbose)
		PRINT_OUT("Step 2 (Algorithm 5):\n");

	Memcpy(subset0, subset1, n);
	while (iter <= *maxiter) {

		// weighted least squares (on return, dat->wx is overwritten by the
		// QR factorization)
		info = fitwls(dat, est->weight, work->dgels_work, work->lwork,
			est->beta, est->resid);
		if (info)
			return WBACON_ERROR_RANK_DEFICIENT;

		// extract L
		for (int i = 0; i < p; i++)
			for (int j = i; j < p; j++)
				est->L[j + i * p] = dat->wx[i + j * n];

		// compute regression scale (sigma)
		est->sigma = regression_scale(est->resid, est->weight, n, p);

		// compute t[i]'s (est->dist)
		err = compute_ti(dat, work, est, subset0, m, est->dist);
		if (err != WBACON_ERROR_OK)
			return err;

		// t-distr. cutoff value (quantile)
		cutoff = qt(*alpha / (double)(2 * (*m + 1)), *m - p, 0, 0);

		// generate new subset that includes all obs. with t[i] < cutoff
		Memcpy(est->weight, dat->w, n);
		*m = 0;
		for (int i = 0; i < n; i++) {
			if (est->dist[i] < cutoff) {
				subset1[i] = 1;
				(*m)++;
			} else {
				subset1[i] = 0;
				est->weight[i] = 0.0;		// weight = 0 if not in subset
			}
		}

		// check whether the subsets differ
		for (i = 0; i < n; i++)
			if (subset0[i] ^ subset1[i])
				break;
		// if the subsets are identical, we return
		if (i == n) {
			*maxiter = iter;
			return WBACON_ERROR_OK;
		}

		if (*verbose)
			PRINT_OUT("  m = %d\n", *m);

		Memcpy(subset0, subset1, n);
		iter++;
	}

	return WBACON_ERROR_CONVERGENCE_FAILURE;
}

/******************************************************************************\
|* select the smallest m observations of array x[n] into the subset           *|
|*  x       array[n]                                                          *|
|*  index   work array[n]                                                     *|
|*  subset  on return: array[n], 1: element is in the subset, 0: otherwise    *|
|*  m       size of the subset                                                *|
|*  n       dimension                                                         *|
\******************************************************************************/
void select_subset(double *x, int *iarray, int *subset, int *m, int *n)
{
	// sort the x[i]'s in ascending order
	psort_array(x, iarray, *n, *m);

	// select the smallest 0...(m-1) observations into the subset
	for (int i = 0; i < *n; i++)
		subset[i] = 0;
	for (int i = 0; i < *m; i++)
		subset[iarray[i]] = 1;
}

/******************************************************************************\
|* Update the cholesky factor and the 'xty' matrix to account for the changes *|
|* between subset0 and subset1                                                *|
|*  dat      typedef struct regdata                                           *|
|*  work     typedef struct workarray                                         *|
|*  est      typedef struct estimate                                          *|
|*  subset0  subset, 1: obs. is in the subset, 0 otherwise, array[n]          *|
|*  subset1  subset, 1: obs. is in the subset, 0 otherwise, array[n]          *|
|*  m        number of obs. in subset1                                        *|
|*  verbose  toggle: 1: additional information is printed to the console;     *|
|*           0: quiet                                                         *|
\******************************************************************************/
wbacon_error_type update_chol_xty(regdata *dat, workarray *work, estimate *est,
	int *subset0, int *subset1, int *verbose)
{
	int n = dat->n, p = dat->p;
	wbacon_error_type err;

	// make copies of L and xty (to restore the arrays if the updating fails)
	Memcpy(work->work_pp, est->L, p * p);
	Memcpy(work->work_np, est->xty, p);

	// first pass: make updates to L and xty (if required)
	int n_update = 0;
	for (int i = 0; i < n; i++) {
		if (subset1[i] > subset0[i]) {
			for (int j = 0; j < p; j++) {
				work->work_p[j] = dat->x[i + j * n] * sqrt(dat->w[i]);
				est->xty[j] += dat->x[i + j * n] * dat->y[i] * dat->w[i];
			}
			chol_update(est->L, work->work_p, p);
			n_update++;
		}
	}

	// in the second pass, we consider the downdates (which may turn L into a
	// rank deficient matrix)
	int n_downdate = 0;
	for (int i = 0; i < n; i++) {
		if (subset1[i] < subset0[i]) {
			for (int j = 0; j < p; j++) {
				work->work_p[j] = dat->x[i + j * n] * sqrt(dat->w[i]);
				est->xty[j] -= dat->x[i + j * n] * dat->y[i] * dat->w[i];
			}
			err = chol_downdate(est->L, work->work_p, p);
			if (err != WBACON_ERROR_OK) {
				// updating failed: restore the original arrays
				Memcpy(est->L, work->work_pp, p * p);
				Memcpy(est->xty, work->work_np, p);
				if (*verbose)
					PRINT_OUT(" (downdate failed, subset is increased)\n");
				return err;
			}
			n_downdate++;
		}
	}
	if (*verbose)
		PRINT_OUT(" (%d up- and %d downdates)\n", n_update, n_downdate);

	return WBACON_ERROR_OK;
}

/******************************************************************************\
|* Rank-one update of the (lower triangular) Cholesky factor                  *|
|*  L   Cholesky factor (lower triangular), array[p, p]                       *|
|*  u   rank-one update, array[p]                                             *|
|*  p   dimension                                                             *|
|*                                                                            *|
|* Golub, G.H, and Van Loan, C.F. (1996). Matrix Computations, 3rd. ed.,      *|
|* Baltimore: The Johns Hopkins University Press, ch. 12.5                    *|
\******************************************************************************/
void chol_update(double *L, double *u, int p)
{
	double a, b, c, tmp;
	for (int i = 0; i < p - 1; i++) {
		tmp = L[i * (p + 1)];						// element L[i,i]
		a = sqrt(_POWER2(tmp) + _POWER2(u[i]));
		b = a / tmp;
		c = u[i] / tmp;
		L[i * (p + 1)] = a;							// element L[i,i]

		for (int j = i + 1; j < p; j++) {			// off-diagonal elements
			L[p * i + j] += c * u[j];
			L[p * i + j] /= b;
			u[j] = b * u[j] - c * L[p * i + j];
		}
	}
	L[p * p - 1] = sqrt(_POWER2(L[p * p - 1]) + _POWER2(u[p - 1]));
}

/******************************************************************************\
|* Rank-one downdate of the (lower triangular) Cholesky factor                *|
|*  L   Cholesky factor (lower triangular), array[p, p]                       *|
|*  u   rank-one update, array[p]                                             *|
|*  p   dimension                                                             *|
|*                                                                            *|
|* NOTE: downdating can turn a full rank matrix into a rank deficient matrix  *|
|* Golub, G.H, and Van Loan, C.F. (1996). Matrix Computations, 3rd. ed.,      *|
|* Baltimore: The Johns Hopkins University Press, ch. 12.5                    *|
\******************************************************************************/
wbacon_error_type chol_downdate(double *L, double *u, int p)
{
	double a, b, c, tmp;
	for (int i = 0; i < p - 1; i++) {
		tmp = L[i * (p + 1)];						// element L[i,i]
		a = _POWER2(tmp) - _POWER2(u[i]);
		if (a < 0.0)
			return WBACON_ERROR_RANK_DEFICIENT;
		a = sqrt(a);
		b = a / tmp;
		c = u[i] / tmp;
		L[i * (p + 1)] = a;							// element L[i,i]

		for (int j = i + 1; j < p; j++) {			// off-diagonal elements
			L[p * i + j] -= c * u[j];
			L[p * i + j] /= b;
			u[j] = b * u[j] - c * L[p * i + j];
		}
	}

	a = _POWER2(L[p * p - 1]) - _POWER2(u[p - 1]);
	if (a < 0.0)
		return WBACON_ERROR_RANK_DEFICIENT;

	L[p * p - 1] = sqrt(a);
	return WBACON_ERROR_OK;
}

/******************************************************************************\
|* Distance measure t[i] of Billor et al. (2000, Eq. 6)                       *|
|*  dat      typedef struct regdata                                           *|
|*  work     typedef struct workarray                                         *|
|*  est      typedef struct estimate                                          *|
|*  subset   subset, 1: obs. is in the subset, 0 otherwise, array[n]          *|
|*  m        number of obs. in subset                                         *|
|*  tis      on entry: work array[n]; on return: t[i]'s                       *|
\******************************************************************************/
wbacon_error_type compute_ti(regdata *dat, workarray *work, estimate *est,
	int *subset, int *m, double* tis)
{
	// compute the diag. of the 'hat' matrix (work->work_n)
	wbacon_error_type err = hat_matrix(dat, work, est->L, work->work_n);
	if (err != WBACON_ERROR_OK)
		return err;

	// compute t[i]'s
	for (int i = 0; i < dat->n; i++) {
		if (subset[i])
			tis[i] = fabs(est->resid[i]) /
				(est->sigma * sqrt(1.0 - work->work_n[i]));
		else
			tis[i] = fabs(est->resid[i]) /
				(est->sigma * sqrt(1.0 + work->work_n[i]));
	}
	return WBACON_ERROR_OK;
}

/******************************************************************************\
|*      *|
\******************************************************************************/
double regression_scale(double *resid, double *w, int n, int p)
{
	double sigma = 0.0, total_w = 0.0;
	for (int i = 0; i < n; i++) {
		total_w += w[i];
		sigma += w[i] * _POWER2(resid[i]);
	}
	sigma /= total_w - p;
	return sqrt(sigma);
}

/******************************************************************************\
|* Least squares estimate (Cholesky factorization)                            *|
|*  L     Cholesky factor (lower triangular), array[p, p]                     *|
|*  x     design matrix, array[n]                                             *|
|*  xty   X^T * y, array[p]                                                   *|
|*  beta  on return: regression coefficients, array[p]                        *|
|*  n, p  dimension                                                           *|
\******************************************************************************/
void cholesky_reg(double *L, double *x, double *xty, double *beta, int *n,
	int *p)
{
	const int int_one = 1;
	double const d_one = 1.0;

	// solve for 'a' (return in beta) in the triangular system L^T * a = xty
	Memcpy(beta, xty, *p);
	F77_CALL(dtrsm)("L", "L", "N", "N", p, &int_one, &d_one, L, p, beta, p);

	// solve for 'beta' in the triangular system L * beta  = a
	F77_CALL(dtrsm)("L", "L", "T", "N", p, &int_one, &d_one, L, p, beta, p);
}

/******************************************************************************\
|* Diagonal elements of the 'hat' matrix                                      *|
|*  dat      typedef struct regdata                                           *|
|*  work     typedef struct workarray                                         *|
|*  L        Cholesky factor (lower triangular), array[p, p]                  *|
|*  hat      diagonal elements of the 'hat' matrix                            *|
\******************************************************************************/
wbacon_error_type hat_matrix(regdata *dat, workarray *work, double *L,
	double *hat)
{
	int n = dat->n, p = dat->p;
	// invert the cholesky factor L
	int info;
	Memcpy(work->work_pp, L, p * p);
	F77_CALL(dtrtri)("L", "N", &p, work->work_pp, &p, &info);
	if (info != 0)
		return WBACON_ERROR_TRIANG_MAT_SINGULAR;

	// triangular matrix multiplication x * L^{-T}
	const double d_one = 1.0;
	Memcpy(work->work_np, dat->x, n * p);
	F77_CALL(dtrmm)("R", "L", "T", "N", &n, &p, &d_one, work->work_pp, &p,
		work->work_np, &n);

	// diagonal elements of hat matrix = row sums of (x * L^{-T})^2
	double tmp;
	for (int i = 0; i < n; i++) {
		tmp = 0.0;
		for (int j = 0; j < p; j++)
			tmp += _POWER2(work->work_np[i + j * n]);
		hat[i] = tmp * dat->w[i];
	}
	return WBACON_ERROR_OK;
}
#undef _POWER2

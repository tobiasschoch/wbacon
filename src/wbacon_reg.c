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

#include "wbacon_reg.h" 
#define _POWER2(_x) ((_x) * (_x))

/* // data structure */
/* typedef struct wbreg_type { */
/* 	int n; */
/* 	int p; */
/* 	double *x;		// design matrix (raw and weighted) */
/* 	double *wx;		 */
/* 	double *y;		// response vector (raw and weighted) */
/* 	double *wy;		 */
/* 	double *beta;	// regression coefficient */
/* 	double *hat;	// 'hat' matrix  */
/* } wbreg */
/*  */
/* // structure of working arrays */
/* typedef struct workarray_struct { */
/* 	int *iarray; */
/* 	double *work_n; */
/* 	double *work_np; */
/* 	double *work_pp; */
/* 	double *work_2n; */
/* } workarray; */

// declarations of local function
wbacon_error_type update_chol_xty(double*, double*, double*, double*, double*,
	double*, double*, int*, int*, int*, int*);
wbacon_error_type chol_downdate(double*, double*, int*);
wbacon_error_type hat_matrix(double*, double*, double*, double*, double*, int*,
	int*);
wbacon_error_type initial_reg(double*, double*, double*, double*, double*, 
	double*, double*, int*, int*, double*, int*, int*, int*, int*); 
void select_subset(double *x, int *iarray, int *subset, int *m, int *n); 
void cholesky_reg(double*, double*, double*, double*, int*, int*);
void chol_update(double*, double*, int*);
void compute_ti(double*, double*, double*, double*, double*, int*, int*, int*,
	int*);


// ----------------------------------------------------------------------------
void print_pp(double *a, int p);
void print_p(double *a, int p);


void test_hat(double *x, double *L, double *hat, int *n, int *p)
{
	wbacon_error_type err;
	double *work_pp = (double*) Calloc(*p * *p, double);
	double *work_np = (double*) Calloc(*n * *p, double);
	
	err = hat_matrix(L, x, hat, work_pp, work_np, n, p);
	if (err != WBACON_ERROR_OK) 
		PRINT_OUT("Error: %s\n", wbacon_error(err));

	Free(work_pp); Free(work_np);
}

void print_pp(double *a, int p)
{
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < p; j++) {
			PRINT_OUT("%.4f\t", a[i + p * j]);
		}
		PRINT_OUT("\n");
	}
	PRINT_OUT("\n");
}

void print_p(double *a, int p)
{
	for (int i = 0; i < p; i++) {
		PRINT_OUT("%.3f\t", a[i]);
	}
	PRINT_OUT("\n\n");
}




/******************************************************************************\
|* BACON regression estimator                                                 *|
\******************************************************************************/
void wbacon_reg(double *x, double *work_x, double *y, double *work_y,
	double *w, double *resid, double *beta, int *subset0, int *iarray,
	double *dist, int *n, int *p, int *m, int *verbose, int *success)
{
	int *subset1;
	double *hat, *L, *work_pp, *tis, *xty, *work_p;

	wbacon_error_type err;
	*success = 1;	

	work_pp = (double*) Calloc(*p * *p, double);
	work_p = (double*) Calloc(*p, double);

	L = (double*) Calloc(*p * *p, double);
	hat = (double*) Calloc(*n, double);
	tis = (double*) Calloc(*n, double);
	subset1 = (int*) Calloc(*n, int);
	xty = (double*) Calloc(*p, double);

	// STEP 0
	// compute the least squares estimate b of the linear system x*b = y, where
	// x[subset0] is of full rank 
	err = initial_reg(x, work_x, y, work_y, w, resid, beta, subset0, iarray, 
		dist, n, p, m, verbose);
	if (err != WBACON_ERROR_OK) {
		*success = 0;
		PRINT_OUT("Error: x %s (step 0)\n", wbacon_error(err));
		goto clean_up;
	}

	// extract R matrix (as a lower triangular matrix: L) 
	for (int i = 0; i < *p; i++) 
		for (int j = i; j < *p; j++) 
			L[j + i * *p] = work_x[i + j * *n];	


print_pp(L, *p);

	// compute xty
	for (int i = 0; i < *p; i++) 
		for (int j = 0; j < *n; j++) 
			if (subset0[j])
				xty[i] += x[j + i * *n] * y[j];	

PRINT_OUT("xty\n");
print_p(xty, *p);

	// compute the diag. of the hat matrix and the t[i]'s 
	err = hat_matrix(L, x, hat, work_pp, work_x, n, p);
	if (err != WBACON_ERROR_OK) {
		*success = 0;
		PRINT_OUT("Error: %s (step 0)\n", wbacon_error(err));
		goto clean_up;
	}

	compute_ti(x, y, tis, beta, hat, subset0, n, p, m);

	// identify the p+1 obs. with the smallest ti's (initial basic subset)
	*m = *p + 1; 
	select_subset(tis, iarray, subset1, m, n); 


//FIXME: 4*p resp. c*p < n (check at caller)
	// STEP 1
	while (*m < *p * 4) {
		// update cholesky factor and xty matrix (subset0 => subset1)
		err = update_chol_xty(x, y, xty, L, work_pp, work_p, work_y, subset0,
			subset1, n, p);

PRINT_OUT("--------------------------------------------------\n");
PRINT_OUT("m = %d\n", *m);
print_pp(L, *p);

PRINT_OUT("xty\n");
print_p(xty, *p);


		if (err != WBACON_ERROR_OK) {
			// enlarge initial subset until L has full rank 
			for (;;) {
				(*m)++;
				subset1[iarray[*m - 1]] = 1;
				// re-do the updating 
				err = update_chol_xty(x, y, xty, L, work_pp, work_p, work_y, 
					subset0, subset1, n, p);
				if (err == WBACON_ERROR_OK)
					break;
				if (*m == *p * 4) {
					*success = 0;
					PRINT_OUT("Error: %s (chol. update, step 1)\n", 
						wbacon_error(err));
					goto clean_up;
				}
			}
			if (*verbose)
				PRINT_OUT("Step 1: subset increased, now m = %d\n", *m);	
		}

		// regression coefficients (with the updated Cholesky factor)
		cholesky_reg(L, x, xty, beta, n, p);		

PRINT_OUT("beta\n");
print_p(beta, *p);

		// compute the t[i]'s (on return: work_y) 
		err = hat_matrix(L, x, hat, work_pp, work_x, n, p);	
		if (err != WBACON_ERROR_OK) {
			*success = 0;
			PRINT_OUT("Error: %s (step 1)\n", wbacon_error(err));
			goto clean_up;
		}

		compute_ti(x, y, work_y, beta, hat, subset1, n, p, m);
		(*m)++;

		// prepare next while iteration
		Memcpy(subset0, subset1, *n);

		// select the m obs. with the smallest t[i]'s
		select_subset(work_y, iarray, subset1, m, n); 
	}

clean_up:
	Free(work_pp); Free(hat); Free(L); Free(tis); Free(subset1); Free(xty);
	Free(work_p);
}

/******************************************************************************\
|* select the smallest m observations of array a[n] into the subset           *|
|*  a       array[n]                                                          *| 
|*  index   work array[n]                                                     *|
|*  subset  on return: array[n], 1: element is in the subset, 0: otherwise    *|
|*  m       size of the subset                                                *|
|*  n       dimension                                                         *|
\******************************************************************************/
void select_subset(double *a, int *iarray, int *subset, int *m, int *n) 
{
	// sort the a[i]'s in ascending order 
	psort_array(a, iarray, *n, *m);

	// select the smallest 0...(m-1) observations into the subset 
	for (int i = 0; i < *n; i++)
		subset[i] = 0;
	for (int i = 0; i < *m; i++) 
		subset[iarray[i]] = 1;
}

/******************************************************************************\
|* Initial basic subset for BACON regression estimator                        *|
\******************************************************************************/
//FIXME: arguments
wbacon_error_type initial_reg(double *x, double *work_x, double *y, 
	double *work_y, double *w, double *resid, double *beta, int *subset, 
	int *iarray, double *dist, int *n, int *p, int *m, int *verbose)
{
	int info;
	double *work, *w_cpy;
	wbacon_error_type status = WBACON_ERROR_OK;
   
	// 1) determine size of work array (with work_x as a dummy for work)
	int lwork = -1;
	fitwls(x, work_x, y, work_y, w, resid, beta, n, p, work_x, &lwork, &info);

	// 2) allocate memory for work array and a copy of w
	work = (double*) Calloc(lwork, double);
	w_cpy= (double*) Calloc(*n, double);

	// 3) compute regression estimate and residuals
	Memcpy(w_cpy, w, *n);
	   
	for (int i = 0; i < *n; i++)
		if (subset[i] == 0)
			w_cpy[i] = 0.0;	    // w = 0 if obs. not in subset

	fitwls(x, work_x, y, work_y, w_cpy, resid, beta, n, p, work, &lwork, &info);

	// if x is rank deficient, we enlarge the subset
	if (info) {
		status = WBACON_ERROR_RANK_DEFICIENT;
		// sort the dist[i]'s in ascending order 
		psort_array(dist, iarray, *n, *n);

		// add obs. to the subset until x has full rank
		while (*m < *n) {
			(*m)++;
			// select obs. with smallest dist (among those obs. not in the
			// subset); and set its weight to the original weight (i.e., not 0)
			int at = iarray[*m - 1];
			subset[at] = 1;
			w_cpy[at] = w[at];
 
			// re-do regression and check rank
			fitwls(x, work_x, y, work_y, w_cpy, resid, beta, n, p, work, &lwork, 
				&info);
			if (info == 0) {
				status = WBACON_ERROR_OK;
				if (*verbose)
					PRINT_OUT("Step 0: subset increased, now m = %d\n", *m);
				break;
			}
		}
	}
	Free(work); Free(w_cpy);
	return status;
}

/******************************************************************************\
|* Update the cholesky factor and the 'xty' matrix to account for the changes *|
|* between subset0 and subset1                                                *|
|*  x        design matrix, array[n, p]                                       *|
|*  y        response variable, array[n]                                      *|
|*  xty      xty matrix, on return: updated, array[p]                         *|
|*  L        Cholesky factor (low. triang.), on return: updated, array[p, p]  *|
|*  work_pp  work array[p, p]                                                 *|
|*  work_p   work array[p]                                                    *|
|*  work2_p  (second) work array[p]                                           *|
|*  subset0  subset with elements {0, 1}, array[n]                            *|
|*  subset1  subset with elements {0, 1}, array[n]                            *|
|*  n, p     dimensions                                                       *|
\******************************************************************************/
wbacon_error_type update_chol_xty(double *x, double *y, double *xty, double *L,
	double *work_pp, double *work_p, double *work2_p, int *subset0, 
	int *subset1, int *n, int *p)
{
	wbacon_error_type err;

	// make copies of L and xty (to restore the arrays if the updating fails) 
	Memcpy(work_pp, L, *p * *p);
	Memcpy(work2_p, xty, *p);

int no_update = 0;
int no_downdate = 0;

	// first pass: make updates to L and xty (if required) 
	for (int i = 0; i < *n; i++) {
		if (subset1[i] > subset0[i]) {
			for (int j = 0; j < *p; j++) {
				work_p[j] = x[i + j * *n];
				xty[j] += x[i + j * *n] * y[i];
			}
			chol_update(L, work_p, p);

no_update++;
		}
	}

PRINT_OUT("updates: %d\n", no_update);


	// in the second pass, we consider the downdates (which may turn L into a 
	// rank deficient matrix)  
	for (int i = 0; i < *n; i++) {
		if (subset1[i] < subset0[i]) {
			for (int j = 0; j < *p; j++) {
				work_p[j] = x[i + j * *n];
				xty[j] -= x[i + j * *n] * y[i];
			}
			err = chol_downdate(L, work_p, p);
			if (err != WBACON_ERROR_OK) {
				// updating failed: restore the original arrays
				Memcpy(L, work_pp, *p * *p);
				Memcpy(xty, work2_p, *p); 

PRINT_OUT("downdate failed\n");
				return err;
			}

no_downdate++;

		}
	}

PRINT_OUT("downdates: %d\n", no_downdate);

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
void chol_update(double *L, double *u, int *p)
{
	double a, b, c, tmp;
	for (int i = 0; i < *p - 1; i++) {
		tmp = L[i * (*p + 1)];						// element L[i,i]
		a = sqrt(_POWER2(tmp) + _POWER2(u[i]));
		b = a / tmp;
		c = u[i] / tmp;
		L[i * (*p + 1)] = a;						// element L[i,i]

		for (int j = i + 1; j < *p; j++) {			// off-diagonal elements
			L[*p * i + j] += c * u[j];
			L[*p * i + j] /= b;
			u[j] = b * u[j] - c * L[*p * i + j];	
		}
	}
	L[*p * *p - 1] = sqrt(_POWER2(L[*p * *p - 1]) + _POWER2(u[*p - 1]));
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
wbacon_error_type chol_downdate(double *L, double *u, int *p)
{
	double a, b, c, tmp;
	for (int i = 0; i < *p - 1; i++) {
		tmp = L[i * (*p + 1)];						// element L[i,i]
		a = _POWER2(tmp) - _POWER2(u[i]);
		if (a < 0.0)
			return WBACON_ERROR_RANK_DEFICIENT;
		a = sqrt(a);
		b = a / tmp;
		c = u[i] / tmp;
		L[i * (*p + 1)] = a;						// element L[i,i]

		for (int j = i + 1; j < *p; j++) {			// off-diagonal elements
			L[*p * i + j] -= c * u[j];
			L[*p * i + j] /= b;
			u[j] = b * u[j] - c * L[*p * i + j];
		}
	}

	a = _POWER2(L[*p * *p - 1]) - _POWER2(u[*p - 1]);
	if (a < 0.0)
		return WBACON_ERROR_RANK_DEFICIENT;

	L[*p * *p - 1] = sqrt(a);
	return WBACON_ERROR_OK;
}

/******************************************************************************\
|* Distance measure t[i] of Billor et al. (2000, Eq. 6)                       *|
|*  x        design matrix, array[n, p]                                       *|
|*  y        response variable, array[n]                                      *|
|*  tis      on entry: work array[n]; on return: t[i]'s                       *|
|*  beta     regression coefficients, array[p]                                *|
|*  hat      diagonal elements of the hat matrix, array[n]                    *|
|*  subset   indicator whether element is in subset or nor, array[n]          *|
|*  n, p     dimensions                                                       *|
|*  m	     size of subset                                                   *|
\******************************************************************************/
void compute_ti(double *x, double *y, double *tis, double *beta, double *hat,
	int *subset, int *n, int *p, int *m)
{
	// compute the residuals (BLAS::dgemv): y = alpha*A*x + beta*y; residuals
	// are returned in 'tis' 
	const int int_one = 1;
	const double d_neg_one = -1.0, d_one = 1.0;
	Memcpy(tis, y, *n);
	F77_CALL(dgemv)("N", n, p, &d_neg_one, x, n, beta, &int_one, &d_one, tis, 
		&int_one);

	// estimate of residual scale
	double sigma = 0.0;
	for (int i = 0; i < *n; i++)
		sigma += _POWER2(tis[i]);

	sigma /= (double)(*m - *p);
	sigma = sqrt(sigma);

	// compute t[i]'s
	for (int i = 0; i < *n; i++) {
		if (subset[i])
			tis[i] = fabs(tis[i]) / (sigma * sqrt(1.0 - hat[i]));
		else
			tis[i] = fabs(tis[i]) / (sigma * sqrt(1.0 + hat[i]));
	}
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
|* Diagonal elements of the hat matrix                                        *|
|*  L        Cholesky factor (lower triangular), array[p, p]                  *|
|*  x	     design matrix, array[n]                                          *|
|*  hat      on return: diagonal elements of hat matrix, array[n]             *|
|*  work_pp  work array[p, p]                                                 *|
|*  work_np  work array[n, p]                                                 *|
|*  n, p     dimension                                                        *|
\******************************************************************************/
wbacon_error_type hat_matrix(double *L, double *x, double *hat,
	double *work_pp, double *work_np, int *n, int *p)
{
	// invert the cholesky factor L  
	int info;
	Memcpy(work_pp, L, *p * *p);
	F77_CALL(dtrtri)("L", "N", p, work_pp, p, &info);
	if (info != 0)
		return WBACON_ERROR_TRIANG_MAT_SINGULAR;

	// triangular matrix multiplication x * L^{-T}
	const double d_one = 1.0;
	Memcpy(work_np, x, *n * *p);
	F77_CALL(dtrmm)("R", "L", "T", "N", n, p, &d_one, work_pp, p, work_np, n);

	// diagonal elements of hat matrix = row sums of (x * L^{-T})^2
	double tmp;
	for (int i = 0; i < *n; i++) {
		tmp = 0.0;
		for (int j = 0; j < *p; j++)
			tmp += _POWER2(work_np[i + j * *n]);
		hat[i] = tmp;
	}
	return WBACON_ERROR_OK;
}
#undef _POWER2 

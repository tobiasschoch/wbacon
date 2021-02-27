/* Weighted linear regression

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
*/

#include "fitwls.h"

/******************************************************************************\
|* Weighted least squares estimate                                            *|
|*                                                                            *|
|*  dat        typedef struct 'regdata'                                       *|
|*  weight     weight used in weighted regression                             *|
|*  work_dgels work array, array[lwork] used for QR factorization in dgels    *|
|*  lwork      size of array 'work_dgels' (if < 0, then 'dgels' determines    *|
|*             the optimal size and returns it as the functions return value) *|
|*  beta       on return, beta is overwritten with the reg. coefficients      *|
|*  resid      on return, resid is overwritten with the residuals             *|
|* NOTE:                                                                      *|
|*  if not successfull 1 is returned; otherwise 0                             *|
|*  dat must contain slots for dat->wx and dat->wy                            *|
|*  on return, dat->wx is overwritten by the QR factorization as returned by  *|
|*  LAPACK: dgeqrf                                                            *|
\******************************************************************************/
int fitwls(regdata *dat, double *weight, double *work_dgels, int lwork,
	double *beta, double *resid)
{
	const int int_1 = 1;
	int info_dgels = 1, n = dat->n, p = dat->p;

	// STEP 0: determine the optimal size of array 'work' and return
	if (lwork < 0) {
		F77_CALL(dgels)("N", &n, &p, &int_1, dat->x, &n, dat->y, &n, work_dgels,
			&lwork, &info_dgels);
		return (int) work_dgels[0];
	}

	// STEP 1: compute least squares fit

	// pre-multiply the design matrix and the response vector by sqrt(w)
	Memcpy(dat->wy, dat->y, n);
	Memcpy(dat->wx, dat->x, n * p);
	double tmp;
	for (int i = 0; i < n; i++) {
		tmp = sqrt(weight[i]);
		dat->wy[i] = dat->y[i] * tmp;

		for (int j = 0; j < p; j++)
			dat->wx[n * j + i] = dat->x[n * j + i] * tmp;
	}

	// weighted least squares estimate (LAPACK::dgels),
	F77_CALL(dgels)("N", &n, &p, &int_1, dat->wx, &n, dat->wy, &n, work_dgels,
		&lwork, &info_dgels);

	// dgels is not well suited as a rank-revealing procedure; i.e., INFO<0
	// iff a diagonal element of the R matrix is exactly 0. This is not
	// helpful; hence, we check the diagonal elements of R separately and
	// issue and error flag if any(abs(diag(R))) is close to zero
	for (int i = 0; i < p; i++)
		if (fabs(dat->wx[(n + 1) * i]) < sqrt(DBL_EPSILON))
			return 1;

	Memcpy(beta, dat->wy, p);		// extract 'beta'

	// residuals
	const double double_minus1 = -1.0, double_1 = 1.0;
	Memcpy(resid, dat->y, n);
	F77_CALL(dgemv)("N", &n, &p, &double_minus1, dat->x, &n, beta, &int_1,
		&double_1, resid, &int_1);

	return 0;
}

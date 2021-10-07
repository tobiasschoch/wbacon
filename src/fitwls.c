/* Weighted linear regression

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
*/

#include "fitwls.h"

/******************************************************************************\
|* Weighted least squares (incl. residual scale estimator and residuals)      *|
|*                                                                            *|
|*  dat        typedef struct 'regdata'                                       *|
|*  est        typedef struct 'estimate'                                      *|
|*  subset     subset of observations                                         *|
|*  work_dgels work array, array[lwork] used for QR factorization in dgels    *|
|*  lwork      size of array 'work_dgels' (if < 0, then 'dgels' determines    *|
|*             the optimal size and returns it as the functions return value) *|
|* NOTE:                                                                      *|
|*  if not successfull 1 is returned; otherwise 0                             *|
|*  dat must contain slots for dat->wx and dat->wy                            *|
|*  on return, dat->wx is overwritten by the QR factorization as returned by  *|
|*  LAPACK: dgeqrf                                                            *|
\******************************************************************************/
int fitwls(regdata *dat, estimate *est, int* restrict subset,
    double* restrict work_dgels, int lwork)
{
    const int int_1 = 1;
    int info_dgels = 1, n = dat->n, p = dat->p;
    double* restrict wx = dat->wx;
    double* restrict wy = dat->wy;
    double* restrict x = dat->x;
    double* restrict y = dat->y;
    double* restrict weight = dat->w;
    double* restrict weight_sqrt = dat->w_sqrt;
    double* restrict beta = est->beta;
    double* restrict resid = est->resid;
    double* restrict sigma = &est->sigma;

    // STEP 0: determine the optimal size of array 'work' and return
    if (lwork < 0) {
        F77_CALL(dgels)("N", &n, &p, &int_1, x, &n, y, &n, work_dgels, &lwork,
            &info_dgels FCONE);
        return (int) work_dgels[0];
    }

    // STEP 1: compute least squares fit
    // pre-multiply the design matrix and the response vector by sqrt(w)
    double sum_w = 0.0, indicator;
    for (int i = 0; i < n; i++) {
        indicator = (double)subset[i];
        sum_w += weight[i] * indicator;
        wy[i] = weight_sqrt[i] * indicator * y[i];
        wx[i] = weight_sqrt[i] * indicator * x[i];
    }

    #pragma omp parallel for if(n > FITWLS_OMP_MIN_SIZE)
    for (int j = 1; j < p; j++) {
        #pragma omp simd
        for (int i = 0; i < n; i++)
            wx[i + n * j] = (double)subset[i] * weight_sqrt[i] * x[i + n * j];
    }

    // weighted least squares estimate (LAPACK::dgels),
    F77_CALL(dgels)("N", &n, &p, &int_1, wx, &n, wy, &n, work_dgels,
        &lwork, &info_dgels FCONE);

    // dgels is not well suited as a rank-revealing procedure; i.e., INFO<0
    // iff a diagonal element of the R matrix is exactly 0. This is not
    // helpful; hence, we check the diagonal elements of R separately and
    // issue and error flag if any(abs(diag(R))) is close to zero
    for (int i = 0; i < p; i++)
        if (fabs(wx[(n + 1) * i]) < sqrt(DBL_EPSILON))
            return 1;

    // extract regression estimates (beta)
    Memcpy(beta, wy, p);

    // residual scale estimate (sigma, using the output of dgels)
    double ssq = 0.0;
    for (int i = p; i < n; i++)
        ssq += wy[i] * wy[i];

    *sigma = sqrt(ssq / (sum_w - (double)p));

    // residuals
    const double double_minus1 = -1.0, double_1 = 1.0;
    Memcpy(resid, y, n);
    F77_CALL(dgemv)("N", &n, &p, &double_minus1, x, &n, beta, &int_1,
        &double_1, resid, &int_1 FCONE);

    return 0;
}

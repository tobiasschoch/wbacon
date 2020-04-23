/*****************************************************************************\
|* sctmv								     *|
|* ------------------------------------------------------------------------- *|
|* PROJECT  sctmv library						     *|
|* SUBEJCT  basic multivariate statistics functions			     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January 19, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include "sctmv.h"

/* some macros */
#define _POWER2(_x) ((_x) * (_x))

/*****************************************************************************\
|*  weighted mean (vector)						     *|
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    x		  array[n, p]				                     *|
|*    w		  array[n], weights			                     *|
|*    center	  on return: array[p]		      	                     *|
|*    ptrn, ptrp  dimensions						     *|
\*****************************************************************************/
void weightedmean(double *x, double *w, double *center, int *ptrn, int *ptrp)
{
   double sum_w = 0.0;
   for (int i = 0; i < *ptrp; i++){
      center[i] = 0.0;
   }

   for (int i = 0; i < *ptrn; i++){
      sum_w += w[i]; 
      for (int j = 0; j < *ptrp; j++){
	 center[j] += x[*ptrn * j + i] * w[i]; 
      }
   }
   for (int j = 0; j < *ptrp; j++){
      center[j] /= sum_w; 
   }
}

/*****************************************************************************\
|*  weighted covariance/ scatter matrix					     *|
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    x		  array[n, p]				                     *|
|*    w		  array[n], weights			                     *|
|*    center	  array[p]		  	      	                     *|
|*    scatter	  on return: array[p, p]	      	                     *|
|*    ptrn, ptrp  dimensions						     *|
|*                                                                           *|
|*  DEPENDENCIES							     *|
|*    DGEMM    LAPACK/BLAS                                                   *|
\*****************************************************************************/
void weightedscatter(double *x, double *w, double *center, double *scatter, 
		     int *ptrn, int *ptrp)
{
   const double done = 1.0, dzero = 0.0;
   double sum_w = 0.0;
   double *x_cpy;
   x_cpy = (double*) Calloc(*ptrn * *ptrp, double); 
   Memcpy(x_cpy, x, *ptrn * *ptrp); 

   // center the data and multiply by sqrt(w) 
   for (int i = 0; i < *ptrn; i++){
      sum_w += w[i];
 
      for (int j = 0; j < *ptrp; j++){
	 x_cpy[*ptrn * j + i] -= center[j];
	 x_cpy[*ptrn * j + i] *= sqrt(w[i]);
      }
   }

   // cross product: scatter matrix = t(x_cpy) %*% x_cpy; 
   F77_CALL(dgemm)("T",	// TRANSA
      "N",		// TRANSB
      ptrp,		// M (number of rows of A)
      ptrp,		// N (number of columns of B)
      ptrn,		// K (number of columns of A)
      &done,		// ALPHA
      x_cpy,		// A 
      ptrn,		// LDA
      x_cpy,		// B
      ptrn,		// LDB
      &dzero,		// BETA
      scatter,		// C (on return: scatter)
      ptrp);		// LDC

   for (int i = 0; i < (*ptrp * *ptrp); i++){
      scatter[i] /= (sum_w - 1); 
   }

   Free(x_cpy); 
}

/*****************************************************************************\
|*  Mahalanobis distance						     *|
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    x		  array[n, p]				                     *|
|*    center	  array[p]	       		      	                     *|
|*    scatter	  array[p, p]		     	      	                     *|
|*    dist	  on return: array[n]	     	      	                     *|
|*    ptrn, ptrp  dimensions						     *|
|*                                                                           *|
|*  DEPENDENCIES							     *|
|*    DPOTRF   LAPACK/BLAS                                                   *|
|*    DTRSM    LAPACK/BLAS                                                   *|
\*****************************************************************************/
void mahalanobis(double *x, double *center, double *scatter, double *dist, 
		 int *ptrn, int *ptrp)
{
   int info;
   double done = 1.0;
   double *x_cpy, *scatter_cpy;

   // center the data 
   x_cpy = (double*) Calloc(*ptrn * *ptrp, double); 
   Memcpy(x_cpy, x, *ptrn * *ptrp); 
   for (int i = 0; i < *ptrn; i++){
      for (int j = 0; j < *ptrp; j++){
	 x_cpy[*ptrn * j + i] -= center[j];
      }
   }

   // Cholesky decomposition of scatter matrix (stored in scatter_cpy) 
   scatter_cpy = (double*) Calloc(*ptrp * *ptrp, double); 
   Memcpy(scatter_cpy, scatter, *ptrp * *ptrp); 
   F77_CALL(dpotrf)("L",   // UPLO
      ptrp,		   // N (dim of A) 
      scatter_cpy,	   // A (on exit L)
      ptrp,		   // LDA
      &info);		
   if (info != 0) Rprintf("mahalanobis: DPOTRF failed\n");

   // Solve for y in A * y = B by forward subsitution (A = Cholesky factor) 
   F77_CALL(dtrsm)("R",	// SIDE
      "L",		// UPLO
      "T",		// TRANSA
      "N",		// DIAG
      ptrn,		// M (number of rows of B)
      ptrp,		// N (number of columns of B)
      &done,		// ALPHA
      scatter_cpy,	// A	
      ptrp,		// LDA
      x_cpy,		// B (on return: y)
      ptrn);		// LDB
 
   // Euclidean norm of y = squared Mahalanobis distance 
   for (int i = 0; i < *ptrn; i++){
      dist[i] = 0.0;
      for(int j = 0; j < *ptrp; j++){
	 dist[i] += _POWER2(x_cpy[*ptrn * j + i]);
      }
      dist[i] = sqrt(dist[i]);
   }

   Free(x_cpy); Free(scatter_cpy);
}

#undef _POWER2 

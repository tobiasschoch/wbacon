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
|*    n, p	  dimensions						     *|
\*****************************************************************************/
void weightedmean(double *x, double *w, double *center, int *n, int *p)
{
   for (int i = 0; i < *p; i++)
      center[i] = 0.0;

   double sum_w = 0.0;
   for (int i = 0; i < *n; i++) {
      sum_w += w[i]; 
      for (int j = 0; j < *p; j++) {
	 center[j] += x[*n * j + i] * w[i]; 
      }
   }
   for (int j = 0; j < *p; j++)
      center[j] /= sum_w; 
}

/*****************************************************************************\
|*  weighted covariance/ scatter matrix					     *|
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    x		  array[n, p]				                     *|
|*    w		  array[n], weights			                     *|
|*    center	  array[p]		  	      	                     *|
|*    scatter	  on return: array[p, p]	      	                     *|
|*    n, p	  dimensions						     *|
|*                                                                           *|
|*  DEPENDENCIES							     *|
|*    DGEMM    LAPACK/BLAS                                                   *|
\*****************************************************************************/
void weightedscatter(double *x, double *w, double *center, double *scatter, 
		     int *n, int *p)
{
   const double done = 1.0, dzero = 0.0;
   double sum_w = 0.0;
   double *x_cpy;

   x_cpy = (double*) Calloc(*n * *p, double); 
   Memcpy(x_cpy, x, *n * *p); 

   // center the data and multiply by sqrt(w) 
   for (int i = 0; i < *n; i++) {
      sum_w += w[i];
 
      for (int j = 0; j < *p; j++) {
	 x_cpy[*n * j + i] -= center[j];
	 x_cpy[*n * j + i] *= sqrt(w[i]);
      }
   }

   // cross product: scatter matrix = t(x_cpy) %*% x_cpy; 
   F77_CALL(dgemm)("T",	// TRANSA
      "N",		// TRANSB
      p,		// M (number of rows of A)
      p,		// N (number of columns of B)
      n,		// K (number of columns of A)
      &done,		// ALPHA
      x_cpy,		// A 
      n,		// LDA
      x_cpy,		// B
      n,		// LDB
      &dzero,		// BETA
      scatter,		// C (on return: scatter)
      p);		// LDC

   for (int i = 0; i < (*p * *p); i++)
      scatter[i] /= (sum_w - 1); 

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
|*    n, p	  dimensions						     *|
|*                                                                           *|
|*  DEPENDENCIES							     *|
|*    DPOTRF   LAPACK/BLAS                                                   *|
|*    DTRSM    LAPACK/BLAS                                                   *|
\*****************************************************************************/
void mahalanobis(double *x, double *center, double *scatter, double *dist, 
		 int *n, int *p)
{
   int info;
   double done = 1.0;
   double *x_cpy, *scatter_cpy;

   // center the data 
   x_cpy = (double*) Calloc(*n * *p, double); 
   Memcpy(x_cpy, x, *n * *p); 
   for (int i = 0; i < *n; i++) {
      for (int j = 0; j < *p; j++) {
	 x_cpy[*n * j + i] -= center[j];
      }
   }

   // Cholesky decomposition of scatter matrix (stored in scatter_cpy) 
   scatter_cpy = (double*) Calloc(*p * *p, double); 
   Memcpy(scatter_cpy, scatter, *p * *p); 
   F77_CALL(dpotrf)("L",   // UPLO
      p,		   // N (dim of A) 
      scatter_cpy,	   // A (on exit L)
      p,		   // LDA
      &info);		
   if (info != 0) Rprintf("mahalanobis: DPOTRF failed\n");

   // Solve for y in A * y = B by forward subsitution (A = Cholesky factor) 
   F77_CALL(dtrsm)("R",	// SIDE
      "L",		// UPLO
      "T",		// TRANSA
      "N",		// DIAG
      n,		// M (number of rows of B)
      p,		// N (number of columns of B)
      &done,		// ALPHA
      scatter_cpy,	// A	
      p,		// LDA
      x_cpy,		// B (on return: y)
      n);		// LDB
 
   // Euclidean norm of y = squared Mahalanobis distance 
   for (int i = 0; i < *n; i++) {
      dist[i] = 0.0;

      for (int j = 0; j < *p; j++) {
	 dist[i] += _POWER2(x_cpy[*n * j + i]);
      }
      dist[i] = sqrt(dist[i]);
   }

   Free(x_cpy); Free(scatter_cpy);
}

#undef _POWER2 

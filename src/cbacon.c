/*****************************************************************************\
|* bacon								     *|
|* ------------------------------------------------------------------------- *|
|* PROJECT  sct library							     *|
|* SUBEJCT  weighted BACON						     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January 19, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include "cbacon.h" 

// macros
#define _RANK_TOLERANCE 1.0e-8

// prototype of local function
static inline double cutoff(double, int, int, int) 
   __attribute__((always_inline));
static inline void initialsubset(double*, double*, double*, int*, int*, int*, 
   int*) __attribute__((always_inline));

/*****************************************************************************\
|*  Chi-squared cutoff value (p-dimensional, sample size n	  	     *|
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    alpha	  scalar (0 <= alpha <= 1)		                     *|
|*    n		  scalar						     *|
|*    p		  scalar	 					     *|
|*                                                                           *|
|*  DEPENDENCIES							     *|
|*    qchisq <R_ext/Utils.h>                                                 *|
|*                                                                           *|
\*****************************************************************************/
static inline double cutoff(double alpha, int n, int k, int p)
{
   double h, chr, cnp, cnpr, cutoff;

   h = ((double)n + (double)p + 1.0) / 2.0;
   chr = fmax(0.0, (h - (double)k) / (h + (double)k));
   cnp = 1.0 + ((double)p + 1.0) / ((double)n - (double)p) + 
      2.0 / ((double)n - 1.0 - 3.0 * (double)p);
   cnpr = cnp + chr;
   cutoff = cnpr * sqrt(qchisq(alpha / (double)n, p, 0, 0));

   return cutoff;
}

/*****************************************************************************\
|* initialsubset	   		     				     *|
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    x		  array[n, p]				                     *|
|*    w		  on exit: elements not in subset have w = 0, array[n]       *|
|*    dist	  array[n]						     *|
|*    subset	  on exit: initial subset		                     *|
|*    ptrn, ptrp  dimensions		  		                     *|
|*									     *| 
|*  DEPENDENCIES							     *|
|*    rsort_with_index	<R_ext/Utils.h>					     *|
|*    DPOTRF LAPACK/ BLAS						     *| 
|*                                                                           *|
\*****************************************************************************/
static inline void initialsubset(double *x, double *w, double *dist, 
   int *subset, int *ptrn, int *ptrp, int *ptrverbose)
{
   int m, info, rank; 
   int *iarray;
   double *dist_sorted, *x_sorted, *w_sorted, *w_sortedcpy, *center, *scatter;

   center = (double*) Calloc(*ptrp, double);
   scatter = (double*) Calloc(*ptrn * *ptrp, double);

   // sort Mahalanobis distance from small to large
   // ------------------------------------------------------------------------
   // make copies because the arrays are sorted in place
   dist_sorted = (double*) Calloc(*ptrn, double);
   Memcpy(dist_sorted, dist, *ptrn);

   x_sorted = (double*) Calloc(*ptrn * *ptrp, double);
   Memcpy(x_sorted, x, *ptrn * *ptrp);

   w_sorted = (double*) Calloc(*ptrn, double);
   Memcpy(w_sorted, w, *ptrn);

   // generate and populate 'iarray'
   iarray = (int*) Calloc(*ptrn, int);
   for (int i = 0; i < *ptrn; i++){
      iarray[i] = i;
   }

   // sort dist (iarray is sorted along with dist_sorted)
   rsort_with_index(dist_sorted, iarray, *ptrn);

   // use iarray to sort x (x_sorted) and w (w_sorted)
   for (int i = 0; i < *ptrn; i++){
      w_sorted[i] = w[iarray[i]];
      for (int j = 0; j < *ptrp; j++){
	 x_sorted[i + *ptrn * j] = x[iarray[i] + *ptrn * j]; 
      }
   }

   // check if scatter matrix of the first 1:m (sorted) observations has full 
   // rank 
   // ------------------------------------------------------------------------
   m = (int)fmin(4.0 * (double)*ptrp, (double)*ptrn * 0.5); 

   w_sortedcpy = (double*) Calloc(*ptrn, double);
   Memcpy(w_sortedcpy, w_sorted, *ptrn);

   while (m < *ptrn){

      // set weights of observations (m+1):n to zero (note: we use int i = m,
      // because of C's zero indexing
      for (int i = m; i < *ptrn; i++){
	 w_sortedcpy[i] = 0.0;	  
      }
      weightedscatter(x_sorted, w_sortedcpy, center, scatter, ptrn, ptrp); 
      
      // Cholesky decomposition of scatter matrix 
      F77_CALL(dpotrf)("L",   // UPLO
	 ptrp,		      // N (dim of A) 
	 scatter,	      // A (on exit L)
	 ptrp,		      // LDA
	 &info);		
      if (info != 0) Rprintf("Initial subset: DPOTRF failed\n");

      // count diagonal elements of Cholesky factor > tol to determine the rank
      // (rank revealing by Cholesky is cheap, but not numerically stable)
      rank = 0;
      for (int i = 0; i < *ptrp; i++){
	 rank += scatter[i * (*ptrp + 1)] > _RANK_TOLERANCE ? 1 : 0;
      }

      if (rank == *ptrp){
	 break;
      }else if (*ptrverbose) {
	 Rprintf("Scatter matrix is rank defficient: subset is enlarged\n");	 
      }

      // prepare next loop
      m++;
      Memcpy(w_sortedcpy, w_sorted, *ptrn);
   }

   // generate initial subset 
   for (int i = 0; i < m; i++){
      subset[iarray[i]] = 1; 
   }

   // set weight = 0 when not in subset
   for (int i = m; i < *ptrn; i++){
      w[iarray[i]] = 0.0;  
   }

   Free(iarray); Free(x_sorted); Free(w_sorted); Free(w_sortedcpy); 
   Free(center); Free(dist_sorted); Free(scatter);
}

/*****************************************************************************\
|*  weighted BACON							     *|
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    x		  array[n, p]                                                *|
|*    w		  array[n]                                                   *|
|*    center	  on return: array[p]                                        *|
|*    scatter	  on return: array[p, p]                                     *|
|*    dist	  on return: array[n]                                        *|
|*    ptrn, ptrp  dimensions                                                 *|
|*    ptralpha	  scalar                                                     *|
|*    init	  array[n]                                                   *|
|*    ptrcutoff	  on return: scalar                                          *|
|*                                                                           *|
|*    DEPENDENCIES                                                           *|
|*	 weightedmeancov                                                     *|
|*       wquantile							     *|
|*       mahalanobis                                                         *|
|*       cutoff								     *|
|*									     *|
\*****************************************************************************/
void wbacon(double *x, double *w, double *center, double *scatter, double *dist, 
   int *ptrn, int *ptrp, double *ptralpha, int *subset, double *ptrcutoff, 
   int *ptrmaxiter, int *ptrverbose)
{
   int converged, subsetsize = 0, iter = 1;
   int *subset0;
   double dhalf = 0.5;
   double *w_cpy, *x_slice; 

   subset0 = (int*) Calloc(*ptrn, int); 
   w_cpy = (double*) Calloc(*ptrn, double); 

   // finding a suitable center and scatter of subset
   // ------------------------------------------------------------------------    
   // coordinate-wise weighted median (initial location) 
   x_slice = (double*) Calloc(*ptrn, double); 
   for (int j = 0; j < *ptrp; j++){
      Memcpy(x_slice, x + *ptrn * j, *ptrn); 
      wquantile(x_slice, w, ptrn, &dhalf, &center[j]);
   } 
   Free(x_slice); 

   // weighted covariance
   weightedscatter(x, w, center, scatter, ptrn, ptrp); 

   // Mahalanobis distance about initial center and scatter 
   mahalanobis(x, center, scatter, dist, ptrn, ptrp);

   // establish inital subset (with the help of the Mahalanobis distance)
   // ------------------------------------------------------------------------    
   Memcpy(w_cpy, w, *ptrn); 
   initialsubset(x, w_cpy, dist, subset, ptrn, ptrp, ptrverbose);

   for (int i = 0; i < *ptrn; i++){
      subsetsize += subset[i];
   }

   for (;;){

      if (*ptrverbose){
	 double percentage;
	 percentage = 100.0 * (double)subsetsize / (double)*ptrn;
	 Rprintf("Subset %d: n = %d (%.1f%%)\n", iter, subsetsize, percentage);
      }

      // weighted mean on the subset 
      weightedmean(x, w_cpy, center, ptrn, ptrp);

      // weighted covariance on the subset 
      weightedscatter(x, w_cpy, center, scatter, ptrn, ptrp);

      // Mahalanobis distance for all elements of the sample
      mahalanobis(x, center, scatter, dist, ptrn, ptrp);

      // convergence check 
      //-----------------------------------------------------------------------
      converged = 0;
      for (int i = 0; i < *ptrn; i++){
	 converged += subset0[i] ^ subset[i]; // check covergence by XOR
      }
      if (converged == 0){
	 *ptrmaxiter = iter;
	 break;    
      }
      // define subset for next iteration 
      //-----------------------------------------------------------------------
      // chi2 cutoff point
      *ptrcutoff = cutoff(*ptralpha, *ptrn, subsetsize, *ptrp);

      // generate new subset
      Memcpy(subset0, subset, *ptrn); 
      Memcpy(w_cpy, w, *ptrn); 

      subsetsize = 0;
      for (int i = 0; i < *ptrn; i++){
	 if (dist[i] < *ptrcutoff){
	    subset[i] = 1;  
	    subsetsize += 1;
	 }else{
	    w_cpy[i] = 0.0;
	 }
      }

      iter++;
      if (iter > *ptrmaxiter) break;
   }
 
   Free(subset0); Free(w_cpy);
}



/*****************************************************************************\
|* sctmv								     *|
|* ------------------------------------------------------------------------- *|
|* PROJECT  sctmv library						     *|
|* SUBEJCT  header file for basic multivariate statistics functions	     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#ifndef _SCTMV_H 
#define _SCTMV_H

// prototypes for the functions
void weightedmean(double*, double*, double*, int*, int*); 
void weightedscatter(double*, double*, double*, double*, int*, int*); 
void mahalanobis(double*, double*, double*, double*, int*, int*);
#endif 

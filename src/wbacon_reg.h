/******************************************************************************\
|* bacon                                                                      *|
|* -------------------------------------------------------------------------- *|
|* PROJECT  sct library                                                       *|
|* SUBEJCT  header file for weighted BACON                                    *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January 19, 2020           *|
|* LICENSE  GPL >= 2                                                          *|
|* COMMENT  [none]                                                            *|
\******************************************************************************/
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "regdata.h"
#include "wquantile.h" 
#include "partial_sort.h"
#include "fitwls.h" 
#include "wbacon_error.h" 

#ifndef _WBACON_REG_H 
#define _WBACON_REG_H

// macros
#define R_PACKAGE 1					// 1: *.dll/*.so for R; 0: standalone binary 
#define _RANK_TOLERANCE 1.0e-8		// criterion to detect rank deficiency

#if R_PACKAGE
	#define PRINT_OUT(_f, ...) Rprintf((_f), ##__VA_ARGS__)
#else
	#define PRINT_OUT(_f, ...) printf((_f), ##__VA_ARGS__)
#endif

// declarations 
void wbacon_reg(double*, double*, double*, double*, double*, int*, double*,
	int*, int*, int*, int*, int*, int*, double*, int*);
#endif 

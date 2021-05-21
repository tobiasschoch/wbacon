#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "wquantile.h"
#include "partial_sort.h"
#include "wbacon_error.h"

#ifdef _OPENMP
    #include <omp.h>
#endif
#define OMP_MIN_SIZE 100000         // OpenMP enabled when n > OMP_MIN_SIZE

#ifndef _WBACON_H
#define _WBACON_H

// macros
#define R_PACKAGE 1					// 1: *.dll/*.so for R; 0: standalone binary
#define _RANK_TOLERANCE 1.0e-8		// criterion to detect rank deficiency

#if R_PACKAGE
    #define PRINT_OUT(_f, ...) Rprintf((_f), ##__VA_ARGS__)
#else
    #define PRINT_OUT(_f, ...) printf((_f), ##__VA_ARGS__)
#endif

// declarations
void wbacon(double*, double*, double*, double*, double*, int*, int*, double*,
    int*, double*, int*, int*, int*, int*, int*, int*);
#endif

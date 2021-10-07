#ifndef  USE_FC_LEN_T
    #define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
    #define FCONE
#endif

#include <Rmath.h>
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

// variadic arguments in macros are supported by gcc (>=3.0), clang (all
// versions), visual studio (>=2005); the version with ##__VA_ARGS__ (which
// will silently eliminate the trailing comma) is not portable (clang complains)
#if R_PACKAGE
    #define PRINT_OUT(...) Rprintf(__VA_ARGS__)
#else
    #define PRINT_OUT(...) printf(__VA_ARGS__)
#endif

// declarations
void wbacon(double*, double*, double*, double*, double*, int*, int*, double*,
    int*, double*, int*, int*, int*, int*, int*, int*);
#endif

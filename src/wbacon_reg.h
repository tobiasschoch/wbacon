#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "regdata.h"
#include "wquantile.h"
#include "partial_sort.h"
#include "fitwls.h"
#include "wbacon_error.h"
#include "median.h"

#ifdef _OPENMP
    #include <omp.h>
#endif
#define REG_OMP_MIN_SIZE 1000000

#ifndef _WBACON_REG_H
#define _WBACON_REG_H

// macros
#define R_PACKAGE 1             // 1: *.dll/*.so for R; 0: standalone binary
#define _RANK_TOLERANCE 1.0e-8  // criterion to detect rank deficiency

#if R_PACKAGE
    #define PRINT_OUT(_f, ...) Rprintf((_f), ##__VA_ARGS__)
#else
    #define PRINT_OUT(_f, ...) printf((_f), ##__VA_ARGS__)
#endif

// declarations
void wbacon_reg(double*, double*, double*, double*, double*, int*, double*,
    int*, int*, int*, int*, int*, int*, double*, int*, int*, int*);
#endif

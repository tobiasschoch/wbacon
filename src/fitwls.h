#ifndef  USE_FC_LEN_T
    #define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
    #define FCONE
#endif

#include <R.h>
#include <Rmath.h>
#include "regdata.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

#ifndef _FITWLS_H
#define _FITWLS_H

#define FITWLS_OMP_MIN_SIZE 1000000

// prototypes for the functions
int fitwls(regdata*, estimate*, int* restrict, double* restrict, int);
#endif

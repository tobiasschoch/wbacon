#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <omp.h>
#include "regdata.h"

#ifndef _FITWLS_H
#define _FITWLS_H

#define FITWLS_OMP_MIN_SIZE 1000000

// prototypes for the functions
int fitwls(regdata*, estimate*, int* restrict, double* restrict, int);
#endif

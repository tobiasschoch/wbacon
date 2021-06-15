/* (partial) sorting of a double array with index

   Copyright (C) 2020-2021 Tobias Schoch (e-mail: tobias.schoch@gmail.com)

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, a copy is available at
   https://www.gnu.org/licenses/

   Reference:  Bentley, J.L. and D.M. McIlroy (1993). Engineering a
               Sort Function, Software - Practice and Experience 23,
               pp. 1249-1265
*/

#include "partial_sort.h"

# define _large_array 50    // pivotal element determined by ninther

static inline void swap2_indx(double* restrict, int* restrict, int, int)
    __attribute__((always_inline));
static inline int med3(double* restrict, int, int, int)
    __attribute__((always_inline));
static inline int choose_pivot(double* restrict, int, int)
    __attribute__((always_inline));
static inline int is_equal(double, double) __attribute__((always_inline));
static inline void partition_3way_indx(double* restrict, int* restrict, int*,
    int*, int*, int*)
    __attribute__((always_inline));
void partial_sort_with_index(double* restrict, int* restrict, int*, int*, int*);

/******************************************************************************\
|* partially sorts an array with index                                        *|
|*  x      on return: partially sorted array[n]                               *|
|*  index  on return: array[n] that is sorted along with x                    *|
|*  n      dimension                                                          *|
|*  k      all elements x[0..k] are sorted into their final position; the     *|
|*         elements (k+1)..n are not sorted                                   *|
\******************************************************************************/
void psort_array(double *x, int *index, int n, int k)
{
    for (int i = 0; i < n; i++)
        index[i] = i;
    int zero = 0, n_minus_one = n - 1, k_minus_one = k - 1;
    partial_sort_with_index(x, index, &zero, &n_minus_one, &k_minus_one);
}

/******************************************************************************\
|* partially sort a double array (along with the index) up to k-th largest    *|
|* element                                                                    *|
|*  x       on return: array[lo..hi] is sorted for the elements x[lo..k]      *|
|*  index   on return: array[lo..hi] is sorted along with x                   *|
|*  lo      dimension (usually: 0)                                            *|
|*  hi      dimension (usually: n - 1)                                        *|
|*  k       integer in 0..(n - 1); determines the k-th largest element up to  *|
|*          which x is to be sorted                                           *|
\******************************************************************************/
void partial_sort_with_index(double* restrict x, int* restrict index,
    int *lo, int *hi, int *k)
{
    if (*hi <= *lo)     // case: n = 1
        return;

//FIXME: insertion sort for small arrays

    // Bentley-McIlroy's 3-way partitioning: the positions of the sentinels 'i'
    // and 'j' are returned
    int i, j;
    partition_3way_indx(x, index, lo, hi, &i, &j);

    // sort all elements that are smaller than the k-th value
    partial_sort_with_index(x, index, lo, &j, k);
    if (*k >= i)
        partial_sort_with_index(x, index, &i, hi, k);
}

/******************************************************************************\
|* Bentley and McIlroy's (1993) 3-way partitioning scheme with index          *|
|*  array   array[lo..hi]                                                     *|
|*  index   array[lo..hi]                                                     *|
|*  lo, hi  dimensions                                                        *|
|*  i, j    sentinels scanning up and down                                    *|
\******************************************************************************/
static inline void partition_3way_indx(double* restrict array,
    int* restrict index, int *lo, int *hi, int *i, int *j)
{
    // determine pivot and swap it into position 'lo' (i.e., position 0)
    swap2_indx(array, index, choose_pivot(array, *lo, *hi), *lo);
    double pivot = array[*lo];

    // Bentley-McIlroy's 3-way partitioning with sentinels i and j,
    // respectively, scanning up and down until they cross; elements equal to
    // the pivot are swapped to the far left and right,

    int p = *lo, q = *hi + 1;
    *i = *lo; *j = *hi + 1;             // initialize the sentinels

    for (;;) {
        while (array[++(*i)] < pivot)
        if (*i == *hi)
            break;

        while (pivot < array[--(*j)])

        if (*j == *lo)
            break;

        if (*i == *j && is_equal(array[*i], pivot))
        swap2_indx(array, index, ++p, *i);

        if (*i >= *j)                   // check if sentinels cross
            break;

        swap2_indx(array, index, *i, *j);

        // swap equal elements to the far left and right, respectively
        if (is_equal(array[*i], pivot))
            swap2_indx(array, index, ++p, *i);

        if (is_equal(array[*j], pivot))
            swap2_indx(array, index, --q, *j);
    }

    // swap equal elements from the borders back to the center
    *i = *j + 1;
    for (int k = *lo; k <= p; k++)
        swap2_indx(array, index, k, (*j)--);

    for (int k = *hi; k >= q; k--)
        swap2_indx(array, index, k, (*i)++);
}

/******************************************************************************\
|* choose pivotal element: for arrays of size < _large_array, the median of   *|
|* three is taken as pivotal element, otherwise we take Tukey's ninther       *|
|*  array   array[lo..hi]                                                     *|
|*  lo, hi  dimension                                                         *|
\******************************************************************************/
static inline int choose_pivot(double* restrict array, int lo, int hi)
{
    int n = hi - lo + 1;
    int mid = lo + n / 2;       // small array: median of three
    if (n > _large_array) {     // large array: Tukey's ninther
        int eps = n / 8;
        lo = med3(array, lo, lo + eps, lo + eps + eps);
        mid = med3(array, mid - eps, mid, mid + eps);
        hi = med3(array, hi - eps - eps, hi - eps, hi);
    }
    return med3(array, lo, mid, hi);
}

/******************************************************************************\
|* swap the elements i and j in array and index                               *|
|*  array  array[n]                                                           *|
|*  index  array[n]                                                           *|
|*  i, j   elements to be swapped                                             *|
\******************************************************************************/
static inline void swap2_indx(double* restrict array, int* restrict index,
    int i, int j)
{
    double tmp0 = array[i]; array[i] = array[j]; array[j] = tmp0;
    int tmp1 = index[i]; index[i] = index[j]; index[j] = tmp1;
}

/******************************************************************************\
|* is_equal: check whether two doubles are equal using Knuth's notion of      *|
|* essential equality                                                         *|
\******************************************************************************/
static inline int is_equal(double a, double b)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) :
        fabs(a)) * DBL_EPSILON);
}

/******************************************************************************\
|* median-of-three (without swaps)                                            *|
\******************************************************************************/
static inline int med3(double* restrict array, int i, int j, int k)
{
    return array[i] < array[j] ?
        (array[j] < array[k] ? j : array[i] < array[k] ? k : i)
        :   (array[j] > array[k] ? j : array[i] > array[k] ? k : i);
}

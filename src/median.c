/* median

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

# define _n_quickselect 40  // switch from insertion sort to quickselect
# define _n_nither 50       // pivotal element determined by ninther

#include "median.h"

static void partition_3way(double *array, int lo, int hi, int *i, int *j);
static inline int choose_pivot(double *array, int lo, int hi)
    __attribute__((always_inline));
static inline void swap(double *array, int i, int j)
    __attribute__((always_inline));
static inline int is_equal(double a, double b) __attribute__((always_inline));
static inline int med3(double *array, int i, int j, int k)
    __attribute__((always_inline));

/******************************************************************************\
|* median (no memory allocation)                                              *|
|*                                                                            *|
|*  array    array[n]                                                         *|
|*  n        dimension                                                        *|
|*  result   on return: median                                                *|
|* NOTE: median destroys the sorting order of array                           *|
\******************************************************************************/
void median_destructive(double *array, int *n, double *result)
{
    int lo = 0, hi = *n - 1;
    int is_even = *n % 2 == 0;
    int k = (*n + 1) / 2 - 1;

    if (*n <= _n_quickselect) {         // insertion sort
        int exch = 0;
        for (int i = hi; i > lo; i--) { // smallest element as sentinel
            if (array[i] < array[i - 1]) {
                swap(array, i, i - 1);
                exch++;
            }
        }
        if (exch != 0) {                // insertion sort with half-exchanges
            int j;
            double pivot;
            for (int i = lo + 2; i <= hi; i++) {
                pivot = array[i];
                j = i;
                while (pivot < array[j - 1]) {
                    array[j] = array[j - 1];
                    j--;
                }
                array[j] = pivot;
            }
        }
    } else {                            // quickselect
        select_k(array, lo, hi, k);
        if (is_even)
            select_k(array, lo, hi, k + 1);
    }

    if (is_even)
        *result = (array[k] + array[k + 1]) / 2.0;
    else
        *result = array[k];
}

/******************************************************************************\
|* select the k-th largest element (for internal use)                         *|
|*                                                                            *|
|*  array   array[lo..hi]                                                     *|
|*  lo      dimension (usually: 0)                                            *|
|*  hi      dimension (usually: n - 1)                                        *|
|*  k       integer in 0:(n - 1)                                              *|
|*                                                                            *|
|* NOTE     select_k sorts 'array' partially, such that element 'k' is in its *|
|*          final (sorted) position => array[k] gives the k-th largest element*|
\******************************************************************************/
void select_k(double *array, int lo, int hi, int k)
{
    if (hi <= lo)       // case: n = 1
        return;

    // Bentley-McIlroy's 3-way partitioning: the positions of the sentinels
    // 'i' and 'j' are returned
    int i, j;
    partition_3way(array, lo, hi, &i, &j);

    // recursion only on the partitioning where element 'k' lies
    if (k <= j)
        select_k(array, lo, j, k);
    else if (k >= i)
        select_k(array, i, hi, k);
}

/******************************************************************************\
|* Bentley and McIlroy's (1993) 3-way partitioning scheme                     *|
|*  array   array[lo..hi]                                                     *|
|*  lo, hi  dimensions                                                        *|
|*  i, j    sentinels scanning up and down                                    *|
\******************************************************************************/
static void partition_3way(double *array, int lo, int hi, int *i, int *j)
{
    // determine pivot and swap it into position 'lo' (i.e., position 0)
    swap(array, choose_pivot(array, lo, hi), lo);
    double pivot = array[lo];

    // Bentley-McIlroy's 3-way partitioning with sentinels i and j,
    // respectively, scanning up and down until they cross; elements equal to
    // the pivot are swapped to the far left and right

    int p = lo, q = hi + 1;
    *i = lo; *j = hi + 1;               // initialize the sentinels

    for (;;) {
        while (array[++(*i)] < pivot)
        if (*i == hi)
            break;

        while (pivot < array[--(*j)])

        if (*j == lo)
            break;

        if (*i == *j && is_equal(array[*i], pivot))
            swap(array, ++p, *i);

        if (*i >= *j)                   // check if sentinels cross
            break;

        swap(array, *i, *j);

        // swap equal elements to the far left and right, respectively
        if (is_equal(array[*i], pivot))
            swap(array, ++p, *i);

        if (is_equal(array[*j], pivot))
            swap(array, --q, *j);
    }

    // swap equal elements from the borders back to the center
    *i = *j + 1;
    for (int k = lo; k <= p; k++)
        swap(array, k, (*j)--);

    for (int k = hi; k >= q; k--)
        swap(array, k, (*i)++);
}

/******************************************************************************\
|* choose pivotal element: for arrays of size < _n_nither, the median of      *|
|* three is taken as pivotal element, otherwise we take Tukey's ninther       *|
|*  array   array[lo..hi]                                                     *|
|*  lo, hi  dimension                                                         *|
\******************************************************************************/
static inline int choose_pivot(double *array, int lo, int hi)
{
    int n = hi - lo + 1;
    int mid = lo + n / 2;       // small array: median of three
    if (n > _n_nither) {        // large array: Tukey's ninther
        int eps = n / 8;
        lo = med3(array, lo, lo + eps, lo + eps + eps);
        mid = med3(array, mid - eps, mid, mid + eps);
        hi = med3(array, hi - eps - eps, hi - eps, hi);
    }
    return med3(array, lo, mid, hi);
}

/******************************************************************************\
|* swap the elements i and j in array                                         *|
|*  array   array[n]                                                          *|
|*  i, j    elements to be swapped                                            *|
\******************************************************************************/
static inline void swap(double *array, int i, int j)
{
    double tmp = array[i];
    array[i] = array[j];
    array[j] = tmp;
}

/******************************************************************************\
|* check whether two doubles are equal using Knuth's notion of essential      *|
|* equality                                                                   *|
\******************************************************************************/
static inline int is_equal(double a, double b)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) :
        fabs(a)) * DBL_EPSILON);
}

/******************************************************************************\
|* median-of-three (but without swaps)                                        *|
\******************************************************************************/
static inline int med3(double *array, int i, int j, int k)
{
    return array[i] < array[j] ?
        (array[j] < array[k] ? j : array[i] < array[k] ? k : i)
        : (array[j] > array[k] ? j : array[i] > array[k] ? k : i);
}

/******************************************************************************\
|* fitwls                                                                     *|
|* -------------------------------------------------------------------------- *|
|* PROJECT  robsurvey library                                                 *|
|* SUBEJCT  header file for utility functions                                 *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), Dec 6, 2020                *|
|* LICENSE  GPL >= 2                                                          *|
|* COMMENT  [none]                                                            *|
\******************************************************************************/
#include <R.h>
#ifndef _UTILS_H
#define _UTILS_H

void print_magic_number(int *set, unsigned int n)
{
    unsigned int at, magic;
    unsigned int count = 0;
    unsigned int loops = n / 10;
    if (n % 10 != 0)
        loops += 1;

    for (int j = 0; j < loops; j++) {
        at = 1;
        magic = 0;
        for (int i = 0; i < 10; i++) {
            magic += at * set[i + j * 10];
            count++;
            if (count == n)
                break;
            at <<= 1;
        }
        Rprintf("%d ", magic);
    }
    Rprintf("\n");
}

#endif

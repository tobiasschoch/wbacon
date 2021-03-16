#include <R.h>
#ifndef _UTILS_H
#define _UTILS_H

/******************************************************************************\
|* prints a 'magic number' that represents an 0-1-array: the number consists  *|
|* of a series/ chunks of base-10 numbers; these chunks represent 10 bits of  *|
|* the 0-1-array (little endian); e.g. '1023' represents the bit pattern/ set *|
|* '1111111111'                                                               *|
|*  set   array[n]                                                            *|
|*  n     dimension                                                           *|
\******************************************************************************/
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

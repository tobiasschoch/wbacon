/* definition of error types

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
*/

#include "wbacon_error.h"

// human readable errors
const char* const CBACON_ERROR_STRINGS[] = {
    "no errors",
    "matrix is rank deficient",
    "matrix is not positive definite",
    "triangular matrix is singular",
    "failure of convergence"
};

// obtain a human readable error message
const char* wbacon_error(wbacon_error_type err)
{
    if (err >= WBACON_ERROR_COUNT)
        return NULL;
    else
        return CBACON_ERROR_STRINGS[err];
}

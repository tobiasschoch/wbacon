#include <stddef.h>

#ifndef _WBACON_ERROR_H
#define _WBACON_ERROR_H

// error handling: error codes
typedef enum wbacon_error_enum {
    WBACON_ERROR_OK = 0,
    WBACON_ERROR_RANK_DEFICIENT,
    WBACON_ERROR_NOT_POSITIVE_DEFINITE,
    WBACON_ERROR_TRIANG_MAT_SINGULAR,
    WBACON_ERROR_CONVERGENCE_FAILURE,
    WBACON_ERROR_COUNT,                 // [not an actual error type]
} wbacon_error_type;

// declaration
const char* wbacon_error(wbacon_error_type);
#endif

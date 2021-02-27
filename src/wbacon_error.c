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

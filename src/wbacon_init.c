// RegisteringDynamic Symbols
#include <R_ext/Rdynload.h>
#include "wbacon.h"
#include "wbacon_reg.h"

// create arrays describing each C routine
static const R_CMethodDef cMethods[]  = {
    {"wbacon", (DL_FUNC) &wbacon, 15},
    {"wbacon_reg", (DL_FUNC) &wbacon_reg, 17},
    {"wquantile", (DL_FUNC) &wquantile, 5},
    {NULL, NULL, 0}
};

// register the C routines to R
void R_init_robsurvey(DllInfo* info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}

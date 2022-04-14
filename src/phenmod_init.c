#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "pimModel.h"

static const R_CMethodDef CEntries[] = {
    {"pimModel", (DL_FUNC) &pimModel, 7},
    {"tsModel",  (DL_FUNC) &tsModel,  5},
    {NULL, NULL, 0}
};

void R_init_phenmod(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


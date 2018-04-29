#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP pcIRT_gamfunk(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"pcIRT_gamfunk", (DL_FUNC) &pcIRT_gamfunk, 2},
    {NULL, NULL, 0}
};

void R_init_pcIRT(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

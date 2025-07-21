#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// 這裡根據你 UpdateBeta 的參數來宣告
void UpdateBeta(double *W, double *Z, double *B, int *dims,
                double *C, double *V, double *alpha,
                double *updateB, int *iteration);

static const R_CMethodDef cMethods[] = {
  {"UpdateBeta", (DL_FUNC) &UpdateBeta, 9},  // 9 是參數個數
  {NULL, NULL, 0}
};

void R_init_PAGE(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

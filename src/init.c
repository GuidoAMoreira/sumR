#include "sumR.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_sumR(DllInfo *info)
{
  /* used by external packages linking to internal sumR code from C */
  R_RegisterCCallable("sumR", "infiniteSum_", (DL_FUNC) &infiniteSum_);
  R_RegisterCCallable("sumR", "infiniteSumToThreshold_",
                      (DL_FUNC) &infiniteSumToThreshold_);
  R_RegisterCCallable("sumR", "infiniteAdaptive_",
                      (DL_FUNC) &infiniteAdaptive_);
  R_RegisterCCallable("sumR", "infiniteCFolding_",
                      (DL_FUNC) &infiniteCFolding_);
}

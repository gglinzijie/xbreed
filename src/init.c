
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/
/* .Fortran calls */
extern void F77_NAME(cf)(int *r_size, int *c_size, int *loci_mat, int *nloci, double *freq1_main);
extern void F77_NAME(hp)(int *hpsize, int *nmarker, int *nqtl, int *loc_c, int *ng, int *in_r_mu, int *in_c_mu, int *in_pois, int *hp_loci, int *method_laf, double *arglaf_F, int *rcs, int *cross);
extern void F77_NAME(ld)(int *r_size, int *c_size, int *loci_mat, int *npair, int *method, double *ld_data);
extern void F77_NAME(sh)(int *in_r_mu, int *in_c_mu, int *in_pois, int *hp_loci, int *in_pois_f, int *hp_loci_f, int *r_rec, int *c_rec, int *nchr, int *szch, int *arg_rec_m, int *r_f, int *c_f, int *rec_f);
static const R_FortranMethodDef FortranEntries[] = {
  {"cf", (DL_FUNC) &F77_NAME(cf),  5},
  {"hp", (DL_FUNC) &F77_NAME(hp), 13},
  {"ld", (DL_FUNC) &F77_NAME(ld),  6},
  {"sh", (DL_FUNC) &F77_NAME(sh), 14},
  {NULL, NULL, 0}
};

void R_init_xbreed(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

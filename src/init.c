// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* void R_init_markovchain(DllInfo* info) { */
/*   R_registerRoutines(info, NULL, NULL, NULL, NULL); */
/*   R_useDynamicSymbols(info, TRUE); */
/* } */

/* .C calls */
extern void coxGL(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void glCoxRate(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void glHaz(void *, void *, void *, void *, void *, void *);
extern void glRate(void *, void *, void *, void *, void *, void *, void *, void *);
extern void glU2(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_ns_est(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _reReg_am1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_re2(SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_reGehan(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_reGehan_s(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_reLog(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_reRate(SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temGehan(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temHaz(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temLog(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temScGehan(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temScLog(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"coxGL",      (DL_FUNC) &coxGL,      11},
    {"glCoxRate",  (DL_FUNC) &glCoxRate,  11},
    {"glHaz",      (DL_FUNC) &glHaz,       6},
    {"glRate",     (DL_FUNC) &glRate,      8},
    {"glU2",       (DL_FUNC) &glU2,        9},
    {"log_ns_est", (DL_FUNC) &log_ns_est, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_reReg_am1",        (DL_FUNC) &_reReg_am1,        6},
    {"_reReg_re2",        (DL_FUNC) &_reReg_re2,        4},
    {"_reReg_reGehan",    (DL_FUNC) &_reReg_reGehan,    5},
    {"_reReg_reGehan_s",    (DL_FUNC) &_reReg_reGehan_s,    6},
    {"_reReg_reLog",      (DL_FUNC) &_reReg_reLog,      5},
    {"_reReg_reRate",     (DL_FUNC) &_reReg_reRate,     4},
    {"_reReg_temGehan",   (DL_FUNC) &_reReg_temGehan,   7},
    {"_reReg_temHaz",     (DL_FUNC) &_reReg_temHaz,     8},
    {"_reReg_temLog",     (DL_FUNC) &_reReg_temLog,     7},
    {"_reReg_temScGehan", (DL_FUNC) &_reReg_temScGehan, 7},
    {"_reReg_temScLog",   (DL_FUNC) &_reReg_temScLog,   7},
    {NULL, NULL, 0}
};

void R_init_reReg(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

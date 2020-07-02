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
extern void alphaEqC(void *, void *, void *, void *, void *, void *);
extern void betaEst(void *, void *, void *, void *, void *, void *, void *, void *);
extern void coxGL(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void glCoxRate(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void glHaz(void *, void *, void *, void *, void *, void *);
extern void glRate(void *, void *, void *, void *, void *, void *, void *, void *);
extern void glU2(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void HWb(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hwHaz(void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_ns_est(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void plLambda(void *, void *, void *, void *, void *, void *, void *);
extern void sarm1(void *, void *, void *, void *, void *, void *, void *);
extern void sc1Gehan(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sc1Log(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sc2(void *, void *, void *, void *, void *, void *, void *);
extern void scRate(void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _reReg_am1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_re2(SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_reGehan(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_reLog(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_reRate(SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temGehan(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temHaz(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temLog(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temScGehan(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _reReg_temScLog(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"alphaEqC",   (DL_FUNC) &alphaEqC,    6},
    {"betaEst",    (DL_FUNC) &betaEst,     8},
    {"coxGL",      (DL_FUNC) &coxGL,      11},
    {"glCoxRate",  (DL_FUNC) &glCoxRate,  11},
    {"glHaz",      (DL_FUNC) &glHaz,       6},
    {"glRate",     (DL_FUNC) &glRate,      8},
    {"glU2",       (DL_FUNC) &glU2,        9},
    {"HWb",        (DL_FUNC) &HWb,         9},
    {"hwHaz",      (DL_FUNC) &hwHaz,       8},
    {"log_ns_est", (DL_FUNC) &log_ns_est, 11},
    {"plLambda",   (DL_FUNC) &plLambda,    7},
    {"sarm1",      (DL_FUNC) &sarm1,       7},
    {"sc1Gehan",   (DL_FUNC) &sc1Gehan,    9},
    {"sc1Log",     (DL_FUNC) &sc1Log,      9},
    {"sc2",        (DL_FUNC) &sc2,         7},
    {"scRate",     (DL_FUNC) &scRate,      9},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_reReg_am1",        (DL_FUNC) &_reReg_am1,        6},
    {"_reReg_re2",        (DL_FUNC) &_reReg_re2,        4},
    {"_reReg_reGehan",    (DL_FUNC) &_reReg_reGehan,    5},
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

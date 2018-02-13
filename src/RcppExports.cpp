// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// trial
DataFrame trial(Rcpp::List categories, std::vector<std::string> genes, std::vector<std::string> pool, int selectSize, int steps, unsigned int randomSeed);
RcppExport SEXP _MulEA_trial(SEXP categoriesSEXP, SEXP genesSEXP, SEXP poolSEXP, SEXP selectSizeSEXP, SEXP stepsSEXP, SEXP randomSeedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type categories(categoriesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type pool(poolSEXP);
    Rcpp::traits::input_parameter< int >::type selectSize(selectSizeSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type randomSeed(randomSeedSEXP);
    rcpp_result_gen = Rcpp::wrap(trial(categories, genes, pool, selectSize, steps, randomSeed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MulEA_trial", (DL_FUNC) &_MulEA_trial, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_MulEA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

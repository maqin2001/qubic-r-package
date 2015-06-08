// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// nothing
NumericMatrix nothing(NumericMatrix matrix);
RcppExport SEXP Qubic_nothing(SEXP matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    __result = Rcpp::wrap(nothing(matrix));
    return __result;
END_RCPP
}
// qubic
List qubic(NumericMatrix matrix);
RcppExport SEXP Qubic_qubic(SEXP matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    __result = Rcpp::wrap(qubic(matrix));
    return __result;
END_RCPP
}

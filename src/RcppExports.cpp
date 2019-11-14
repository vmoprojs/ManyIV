#include <RcppArmadillo.h>
#include <RcppEigen.h>

/*https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
 */

 using namespace Rcpp;


SEXP arMM(  arma::mat A,   arma::mat B);
RcppExport SEXP _ManyIV_arMM(SEXP nA, SEXP nB) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter<  arma::mat >::type A(nA);
    Rcpp::traits::input_parameter<  arma::mat >::type B(nB);
    rcpp_result_gen = Rcpp::wrap(arMM(A,B));
    return rcpp_result_gen;
END_RCPP
}



SEXP eiMM(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _ManyIV_eiMM(SEXP nA, SEXP nB) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter<  Eigen::MatrixXd >::type A(nA);
    Rcpp::traits::input_parameter<  Eigen::MatrixXd >::type B(nB);
    rcpp_result_gen = Rcpp::wrap(eiMM(A,B));
    return rcpp_result_gen;
END_RCPP
}


SEXP eiMu(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _ManyIV_eiMMM(SEXP nA, SEXP nB) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(nA);
    Rcpp::traits::input_parameter<  Eigen::Map<Eigen::MatrixXd> >::type B(nB);
    rcpp_result_gen = Rcpp::wrap(eiMu(A,B));
    return rcpp_result_gen;
END_RCPP
}

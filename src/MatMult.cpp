#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "header.h"



// [[Rcpp::export]]
SEXP arMM( arma::mat A,  arma::mat B){
  arma::mat C = A * B;
  return Rcpp::wrap(C);
}



// [[Rcpp::export]]
SEXP eiMM(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eiMu(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}




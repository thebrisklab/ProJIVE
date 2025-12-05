// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Eigen::MatrixXd demo_inverse(const Eigen::Map<Eigen::MatrixXd> & x) {

  Eigen::MatrixXd x_inverse = x.inverse();

  return x_inverse;

}



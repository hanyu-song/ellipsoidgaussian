// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>


//' Matrix multiplication
//'
//' @description
//' `armaMatMult` calculates the matrix product of A and B.
//'
//' @param A a matrix
//' @param B a matrix
//'
// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
  arma::mat C = A * B;

  return Rcpp::wrap(C);
}

//' Matrix multiplication
//'
//' @description
//' `eigenMatMult` calculates the matrix product of A and B.
//'
//' @param A a matrix
//' @param B a matrix
//'
// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}

//' Matrix multiplication
//'
//' @description
//' `eigeMapnMatMult` calculates the matrix product of A and B.
//'
//' @param A a matrix
//' @param B a matrix
//'
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//[[Rcpp::export]]
arma::mat rotate_cpp(const arma::vec &a, const arma::vec& b) {
  //int p = size(a);
  double ab = sum(a.t()*b);
  arma::mat ca = a - b*ab;
  //ca = NumericMatrix(p, 1, ca.begin());
  double scaler = sqrt(accu(ca%ca));
  ca = ca/scaler;
  arma::mat A = b * ca.t();
  A = A - A.t();
  double theta = acos(ab);
  arma::mat D(size(A), arma::fill::eye);
  arma::mat R = D + sin(theta) * A + (cos(theta) - 1) * (b * b.t() + ca * ca.t());
  return R;
}

#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <utility>
#include <boost/math/tools/roots.hpp>
#include <stdlib.h>
#include <Rcpp.h>
class tolerance {
public:
  tolerance(double eps) :
  _eps(eps) {
  }
  bool operator()(double a, double b) {
    return (fabs(b - a) <= _eps);
  }
private:
  double _eps;
};


//' Helper function for evaluating the pseudo-normalising constant
//'
//' @description
//' `kfb_cpp` is a helper function for evaluating the pseudo-normalising constant.
//' See Kume and Wood (2005) for more details.
//'
//' @references
//' Kume, A., and Andrew T. A. Wood. “Saddlepoint Approximations for the Bingham
//' and Fisher-Bingham Normalising Constants.” Biometrika 92, no. 2 (2005):
//' 465–76. http://www.jstor.org/stable/20441200.
//'
//' @param j The power
//' @param gam Gamma.
//' @param lam Lambda.
//' @param ta t.
//'
//[[Rcpp::export]]
double kfb_cpp(const unsigned& j, const arma::vec& gam, const arma::vec& lam, double ta) {
  double kd;
  if (j == 1) {
    kd = sum(0.5/(lam - ta) + 0.25 * pow(gam,2) /pow(lam - ta, 2));
  }
  else if (j > 1){
    kd = sum(0.5 * std::tgamma(j)/pow(lam - ta, j) + 0.25 *
      std::tgamma(j+1) * pow(gam, 2)/pow(lam - ta,j + 1));
  }
  else {
    Rcpp::stop("Input j is invalid!\n");
  }
  return(kd);
}


//' Find saddle point
//' @description
//' `saddle_equat_cpp` is a helper function for evaluating the pseudo-normalising constant.
//' See Kume and Wood (2005) for more details.
//'
//' @references
//' Kume, A., and Andrew T. A. Wood. “Saddlepoint Approximations for the Bingham
//' and Fisher-Bingham Normalising Constants.” Biometrika 92, no. 2 (2005):
//' 465–76. http://www.jstor.org/stable/20441200.
//'
//' @param ta t.
//' @param gam Gamma.
//' @param lam Lambda.
//'
//[[Rcpp::export]]
double saddle_equat_cpp(const double& ta, const arma::vec& gam,const arma::vec& lam) {
  arma::uword len = gam.n_elem;
 // assert(len == lam.n_elem); for speed purpose
  arma::vec term2(len);
  double const_ = log(4);
  double res = 0.0;
  arma::vec temp(len);
  temp = lam - ta;
  term2 = 2 * log(abs(gam)) -  2 * log(abs(temp)) - const_;
 // term2 = log(pow(gam, 2)) - log(pow(lam - ta, 2));
 // term2 = log(pow(gam, 2)) - log(pow(temp, 2)) - const_;
  //res = sum(0.5 / (lam - ta) + exp(term2)) - 1.0;
 // res = sum(0.5 / temp + exp(term2)) - 1.0;
  res = 0.5 * sum(1.0 / temp) + sum(exp(term2))  - 1.0;
  return(res);
}

//' Reorder vector
//'
//' @description
//' `reorder_vec` orders vector \code{v} based on the order \code{idx}.
//'
//' @param v A vector.
//' @param idx A vector of indices.
//'
//[[Rcpp::export]]
arma::vec reorder_vec(const arma::vec &v, const arma::uvec &idx) {
  arma::uword len = v.n_elem;
  arma::vec res(len);
  for (arma::uword i = 0; i < len; ++i) {
    res[i]= v[idx[i]];
  }
  return(res);
}

//' Find the root to the saddle point equation.
//'
//' @description
//' `root4SaddleEquat` finds the root to the saddlepoint equation
//'
//' @param gam A vector of length n.
//' @param lam A vector of length n.
//' @param low The lower bound of the solution.
//' @param up The upper bound of the solution.
//'
//[[Rcpp::export]]
double root4SaddleEquat(const arma::vec& gam, const arma::vec& lam,
                        double low, double up) {
  boost::uintmax_t max_iter = 1000;
  double res;
  tolerance tol = 1e-6;
  std::pair<double, double> found =
    boost::math::tools::toms748_solve([gam, lam](double x){ return saddle_equat_cpp(x, gam, lam);},
      low, up, tol,max_iter);
  while (max_iter > 1000) { // did not test this part
    low = found.first;
    up = found.second;
    max_iter = 1000;
    found =
      boost::math::tools::toms748_solve([gam, lam](double x){ return saddle_equat_cpp(x, gam, lam);},
        low, up, tol,max_iter);
  }
  res = found.first;
  return(res);
}


//' Fisher-Bingham normalising constant
//'
//' @description
//' `findFBconst_cpp` calculates the normalising constant in the Fisher-Bingham
//' distribution.
//'
//' @param gam A vector
//' @param lam A vector
//' @param which_ The order of approximation. 1: First order; 2: second order first
//' type; 3: second order second type.
//' @param ordered Whether gam or lam is ordered.
//'
//[[Rcpp::export]]
double findFBconst_cpp(arma::vec& gam, arma::vec& lam, const arma::uword& which_, const bool& ordered) {
  double res;
  if(!ordered) {
    arma::uvec ord = arma::stable_sort_index(lam);
    gam = reorder_vec(gam, ord);
    lam = reorder_vec(lam, ord);
  }
  arma::uword p = gam.n_elem;
  double mina = lam[0];
  double aaa = 0;
  if (mina <= 0) {
    aaa = -mina + 1;
    lam = lam + aaa;
  }
  double gam_max =  max(abs(gam));
  double low = lam[0] - 0.25 * p - 0.5 * sqrt(0.25 * pow(p,2) + p *  pow(gam_max,2));
  double up = lam[0] - 0.25 - 0.5 * sqrt(0.25 + pow(gam[0],2));
  double tau = root4SaddleEquat(gam, lam, low, up);
  double temp = kfb_cpp(2, gam, lam, tau);
  double rho3 = kfb_cpp(3, gam, lam, tau)/pow(temp,1.5);
  double rho4 = kfb_cpp(4, gam, lam, tau)/pow(temp, 2);
  double Ta = rho4 / 8.0 - 5.0 / 24.0 * pow(rho3, 2);
  double c1 = 0.5 * log(2) + 0.5 * (p - 1) * log(arma::datum::pi) -
    0.5 * log(temp) - 0.5 * sum(log(lam - tau)) - tau + 0.25 *
    sum(pow(gam,2)/(lam - tau));
  if (which_ == 1) {
    res = c1;
  }
  else if (which_ == 2) {
    res = c1 + log1p(Ta);
  }
  else {
    res = c1 + Ta;
  }

  if (mina <= 0) {
    res += aaa;
  }
  //arma::vec logcon = {c1, c2, c3};
  return(res);
}



//' Fisher-Bingham normalising constant
//'
//' @description
//' `approxFBconst_cpp` approximates the normalising constant in the Fisher-Bingham
//' distribution using saddlepoint approximation (Kume and Wood 2005)
//'
//' @param para1 The vector parameter, length is k.
//' @param para2 The matrix parameter, k by k.
//' @param idx_ The order of approximation,1: First order; 2: second order first
//' type; 3: second order second type.
//'
//' @references
//'  Kume, A., and Andrew T. A. Wood. “Saddlepoint Approximations for the Bingham
//'  and Fisher-Bingham Normalising Constants.” Biometrika 92, no. 2 (2005):
//'  465–76. http://www.jstor.org/stable/20441200.
//[[Rcpp::export]]
arma::vec approxFBconst_cpp(const arma::mat& para1, const arma::mat& para2,
                     const arma::uword& idx_) {
  arma::uword n = para1.n_cols; // k by n
// std::vector<double> res(n); return type can be:std::vector<double> instead
  arma::vec res(n);
  bool ordered_ = 1;
  arma::vec eigVals;
  arma::mat eigVec;
  arma::eig_sym(eigVals, eigVec, para2);
// arma::uvec ord = arma::stable_sort_index(eigVals);
 //no need to reorder eigVals, since eigVals are in acsending order
 for (arma::uword i = 0; i < n; ++i) {
   arma::vec gam = eigVec.t() * para1.col(i);
   res[i] = findFBconst_cpp(gam, eigVals, idx_, ordered_);
 }
  return(res);
}

//' Convert a vech format to a matrix
//'
//' @description
//' `Vech2Mat_cpp` converts a single vector of length \eqn{d(d+1)/2} obtained
//' through the vech (vector half) operation to its original matrix format. The vech
//' (vector half) operator takes a symmetric \eqn{d \times d} matrix and
//' stacks the lower triangular half into a single vector of length \eqn{d(d+1)/2}.
//'
//' @param para2_vech A vector obtained through vech operation.
//'
//[[Rcpp::export]]
arma::mat Vech2Mat_cpp(const arma::vec& para2_vech) {
  arma::uword n = (-1 + sqrt(1 + 8 * para2_vech.n_elem)) / 2;
  arma::mat V(n,n);
  // make empty matrices
  arma::mat Z(n,n,arma::fill::zeros);
 arma::mat X(n,n,arma::fill::zeros);
  // fill matrices with integers
  arma::vec idx = arma::linspace<arma::mat>(1,n,n);
  X.each_col() += idx;
  Z.each_row() += trans(idx);
  // assign upper triangular elements
  // the >= allows inclusion of diagonal elements
  V.elem(find(Z<=X)) = para2_vech;
 /*** arma::uword idx = 0;
  for (arma::uword j = 0; j < n; ++j) { // col.
    for (arma::uword i = j; i < n; ++i) { // row.
      V(i,j) = para2_vech[idx];
      idx++;
    }
  }
  ****/
  V = arma::symmatl(V);
  return(V);
}


//' Gradient w.r.t. the matrix parameter in the pseudo-normalising constant
//'
//' @description
//' `calclogPseudoconst_MatParGrad4_cpp` calculates the gradient of the logarithm
//' of the pseudo-normalising constant w.r.t. the matrix parameter.
//'
//' @param para2_vech A vector obtained through vech operation from the matrix
//' parameter.
//' @param para1 A vector, the vector parameter in the pseudo-normalising constant.
//'
//[[Rcpp::export]]
double calclogPseudoconst_MatParGrad4_cpp(const arma::vec& para2_vech, const arma::mat& para1) {
  // para1: k by n;
  // p <- length(invSigma); k <- nrow(para1)
  //  lambda2 <- matrix(lambda, nrow = p, ncol = k)
  // para2 <- t(lambda2) %*% diag(invSigma) %*% lambda2 / 2
  // lam: eigenvalue
  // gam: vector
  //#para2 <- rockchalk::vech2mat(para2_vech, lowerOnly = F)
  arma::mat para2 = Vech2Mat_cpp(para2_vech);
  double res = sum(approxFBconst_cpp(para1, para2, 3));
  return(res);
}

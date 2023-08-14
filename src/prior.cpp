// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>


double calcloglaplaceDeri_cpp(const double& var,
                       const double& scale,
                       const double& mean) {
  double res = 1 / scale;
  if (var > mean) {
    res = -res;
  }
  return(res);
}

//' Prior gradient of Lambda
//'
//' @description
//' `calcgradLambdaprior_cpp` calculates the gradient of the Dirichlet-Laplace
//' prior w.r.t. Lambda,  a factor loading matrix.
//'
//' @param Lambda A numeric matrix, the factor loading matrix.
//' @param phi A prior parameter in the Dirichlet-laplace prior.
//' @param kappa A prior parameter in the Dirichlet-Lapalce prior.
//'
//[[Rcpp::export]]
std::vector<double> calcgradLambdaprior_cpp(const arma::mat& Lambda,
                                        const arma::mat& phi,
                                        const arma::vec& kappa) {
  arma::uword p = Lambda.n_rows;
  arma::uword k = Lambda.n_cols;
  arma::vec temp(p * k);
  arma::uword pos = 0;
  for (arma::uword h = 0; h < k; ++h) {
    for (arma::uword j = 0; j < p; ++j) {
      temp[pos] = calcloglaplaceDeri_cpp(Lambda(j,h), phi(j,h) * kappa[j], 0);
      pos++;
    }
  }
  std::vector<double> res = arma::conv_to<std::vector<double>>::from(temp);
    return(res);
}

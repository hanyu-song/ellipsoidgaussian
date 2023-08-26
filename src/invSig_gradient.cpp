// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <Rcpp.h>


//' Lower Triangular of the outer product
//'
//' @description
//' `lowerTriOuterProduct_cpp` calculates the outer product of \code{vec}
//' and returns the lower triangular part of the resulting matrix
//'
//' @param vec A numeric vector.
//'
//[[Rcpp::export]]
arma::rowvec lowerTriOuterProduct_cpp(const arma::rowvec& vec) {
// take the lower triangular of the outer product
  arma::uword k = vec.n_elem;
  arma::rowvec res(k * (k + 1) / 2);
  arma::uword idx = 0;
  for (arma::uword j = 0; j < k; ++j) {
    for (arma::uword i = j; i < k; ++i){
      res[idx] = vec[i] * vec[j];
      idx = idx + 1;
    }
  }
  return(res);
}
//' Gradient change due to log transformation
//'
//' @description
//' `logTransGradAdjust_cpp` calculates the gradient change due to log transformation.
//'
//' @param invSig A numeric vector, the parameter values.
//' @param invSig_jaco The associated Jacobian.
//'
//[[Rcpp::export]]
arma::vec logTransGradAdjust_cpp(const arma::vec& invSig, const arma::vec& invSig_jaco) {
  arma::vec res= invSig % invSig_jaco;
  return(res);
}

//' Gradient for the precision parameters
//'
//' @description
//' `calcPsdMat_lingSiv_jaco_cpp` calculates the Jacobian of the pseudo-normalising constant
//' w.r.t. \code{invSig} through change in the positive semi-definite matrix parameter.
//'
//' @param invSig A numeric vector, the precision parameters.
//' @param psdMat_grad The Jacobian w.r.t. the positive semi-definite matrix.
//' @param lambda A numeric matrix of p by k, the factor loading matrix.
//[[Rcpp::export]]
arma::vec calcPsdMat_linvSig_jaco_cpp(const arma::vec& invSig,
                                      const arma::vec& psdMat_grad,
                                      const arma::mat& lambda) {
// psdMat_grad: k(k+1)/2
  arma::uword p = lambda.n_rows;
  arma::uword k = lambda.n_cols;
  arma::mat jaco(p, k * (k + 1) / 2);
  for(arma::uword i = 0; i < p; ++i) {
   // arma::mat temp = trans(lambda.row(i, )) * lambda.row(i,); // k by k
   // jaco.row(i) = temp.elem(arma::trimatl_ind(arma::size(temp))) / 2;
   jaco.row(i) = lowerTriOuterProduct_cpp(lambda.row(i))/ 2;
  }
  arma::vec invSig_jaco = jaco * psdMat_grad;
  arma::vec res = logTransGradAdjust_cpp(invSig, invSig_jaco);
  return(res);
}


//' Gradient for the precision parameters
//'
//' @description
//' `calcVec_linvSig_jaco_cpp` calculates the Jacobian of the pseudo-normalising constant
//' w.r.t. \code{invSig} through change in the vector parameter.
//'
//' @param invSig A numeric vector the precision parameters.
//' @param vec_grad The gradient of the pseudo-normalising constant w.r.t. the vector.
//' @param lambda A numeric matrix, the factor loading matrix.
//' @param centered_sdat A numeric matrix, centered minibatch.
//'
//[[Rcpp::export]]
arma::vec calcVec_linvSig_jaco_cpp(const arma::vec& invSig,
                                   const arma::mat& vec_grad,
                                   const arma::mat& lambda,
                                   const arma::mat& centered_sdat) {
//vec_grad: k by n
// lambda
arma::uword p = lambda.n_rows;
 // arma::uword k = lambda.n_cols;
 // arma::uword n = centered_sdat.n_rows;
  arma::vec res(p);
  for (arma::uword j = 0; j < p; ++j) {
    //temp = trans(centered_sdat.col(j) * trans(lambda.row(j)));
    arma::mat temp = trans(lambda.row(j)) * trans(centered_sdat.col(j)); // k by n
    res[j] = accu(temp % vec_grad);
  }
  res = logTransGradAdjust_cpp(invSig, res);
  return(res);
}

//' Gradient for the precision parameters
//'
//' @description
//' `calc_linvSig_grad_cpp` calculates the gradient of the log-likelihood w.r.t.
//' the function of the precisions that is being updated in the sampling procedure.
//'
//' @param invSig The precision parameters.
//' @param psdMat_grad The gradient w.r.t. the matrix parameter in the pseudo-normalising
//' constant.
//' @param vec_grad The gradient w.r.t. the vector parameter in the pseudo-normalising
//' constant.
//' @param lambda A numeric matrix, the factor loading matrix.
//' @param centered_sdat A numeric matrix, centered minibatch data.
//'
//[[Rcpp::export]]
std::vector<double> calc_linvSig_grad_cpp(const arma::vec& invSig,
                                      const arma::vec& psdMat_grad,
                                      const arma::mat& vec_grad,
                                      const arma::mat& lambda,
                                      const arma::mat& centered_sdat) {
 arma::vec temp = calcPsdMat_linvSig_jaco_cpp(invSig, psdMat_grad, lambda) +
calcVec_linvSig_jaco_cpp(invSig, vec_grad, lambda, centered_sdat);
  std::vector<double> res = arma::conv_to<std::vector<double>>::from(temp);
  return(res);
}


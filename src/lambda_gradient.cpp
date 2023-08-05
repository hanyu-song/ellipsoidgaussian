// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]] 
#include <RcppArmadillo.h>
#include <cmath>
#include <Rcpp.h>
//[[Rcpp::export]]
arma::mat calc_lambda_grad_MatPart2_cpp(const arma::mat& lambda, 
                                    const arma::vec& psdMat_grad, 
                                    const arma::vec& invSigma) {
  int p = lambda.n_rows;
  arma::uword k = lambda.n_cols;
  arma::mat mat_grad3(k * (k + 1) / 2, p * k);
  for(int j = 0; j < k; ++j) {
    for(int i = 0; i < p; ++i) {
      arma::uword pos = j * p + i;
      arma::mat temp; temp.zeros(k, k);
      temp.col(j) = invSigma[i] * lambda.row(i).t() / 2;
      temp.row(j) = temp.row(j) + invSigma[i] * lambda.row(i) / 2;
      arma::vec temp2 = temp.elem(arma::trimatl_ind(size(temp)));
      mat_grad3.col(pos) = temp2;
    }
  }
 // arma::mat temp3 = psdMat_grad.t() * mat_grad3; // 1 * (pk)
 // std::vector<double> part1 = arma::conv_to<std::vector<double>>::from(temp3);
 arma::mat part1 = psdMat_grad.t() * mat_grad3;;
  return(part1);
}


//[[Rcpp::export]]
arma::mat calc_lambda_grad_VecPart_cpp(const arma::vec& invSigma, 
                                             const arma::mat& centered_sdat,
                                             const arma::mat& vec_grad) {
  arma::uword nsdat = centered_sdat.n_rows;
  arma::uword p = centered_sdat.n_cols;
  arma::uword k = vec_grad.n_rows;
  arma::mat part2;
  part2.zeros(1,p * k);
  std::vector<double> res;
  arma::uword pos;
  for (arma::uword l = 0; l < nsdat; ++l) {
    arma::mat temp;
    temp.zeros(k, p * k); // k by pk
    for (arma::uword j = 0; j < k; ++j) {
      for(arma::uword i = 0; i < p; ++i) {
        pos = j * p + i;
    //    Rcpp::Rcout << pos << '\n';
        temp(j, pos) = invSigma[i] * centered_sdat(l, i);
        }
    }
    arma::mat temp3 = (vec_grad.col(l)).t() * temp; // 1 * (pk)
    part2 = part2 + temp3;
  }
//  res = arma::conv_to<std::vector<double>>::from(part2);
//return(res);
return(part2);
}

//[[Rcpp::export]]
std::vector<double> calc_lambda_grad_cpp(const arma::mat& lambda,
                                         const arma::vec& psdMat_grad, 
                                         const arma::mat& vec_grad,
                                         const arma::vec& invSigma, 
                                         const arma::mat& centered_sdat) {
  arma::mat part2 = calc_lambda_grad_VecPart_cpp(invSigma, centered_sdat, vec_grad);
  arma::mat part1 = calc_lambda_grad_MatPart2_cpp(lambda, psdMat_grad, invSigma);
  arma::mat temp = part1 + part2;
  std::vector<double> res = arma::conv_to<std::vector<double>>::from(temp);
  return(res);
}
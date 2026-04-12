#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::mat;

// [[Rcpp::export]]
SEXP dense_mat_xptr_cpp(const arma::mat& X) {
  return Rcpp::XPtr<arma::mat>(new arma::mat(X), true);
}

// [[Rcpp::export]]
arma::mat dense_matmul_xptr_cpp(SEXP A_ptr, const arma::mat& B) {
  Rcpp::XPtr<arma::mat> ptr(A_ptr);
  const arma::mat& A = *ptr;
  if (A.n_cols != B.n_rows) {
    Rcpp::stop("incompatible matrix dimensions");
  }
  return A * B;
}

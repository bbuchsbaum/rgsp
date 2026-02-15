#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::sp_mat;
using arma::mat;
using arma::vec;
using arma::uword;

// Lanczos approximation for f(L) X where f is provided as an R function.
// K steps (truncated if convergence reached). Handles multiple signal columns.
// [[Rcpp::export]]
arma::mat lanczos_filter_cpp(const arma::sp_mat& L,
                             const arma::mat& X,
                             int K,
                             Rcpp::Function kernel) {
  const uword n = L.n_rows;
  if (L.n_cols != n) Rcpp::stop("L must be square");
  if (X.n_rows != n) Rcpp::stop("X rows must match L dimension");
  if (K < 1) Rcpp::stop("K must be >= 1");

  mat Y(n, X.n_cols, arma::fill::zeros);
  const int Kmax = std::min<int>(K, static_cast<int>(n));

  // Reusable buffers
  vec v(n, arma::fill::zeros);
  vec v_prev(n, arma::fill::zeros);
  vec w(n, arma::fill::zeros);
  mat V(n, Kmax, arma::fill::zeros);
  std::vector<double> alphas;
  std::vector<double> betas;
  alphas.reserve(Kmax);
  betas.reserve(Kmax);

  for (uword col = 0; col < X.n_cols; ++col) {
    v = X.col(col);
    double normx = arma::norm(v, 2);
    if (normx == 0.0) {
      continue; // zero column stays zero
    }

    v /= normx;
    v_prev.zeros();
    alphas.clear();
    betas.clear();

    int keff = 0;
    for (int k = 0; k < Kmax; ++k) {
      V.col(k) = v;
      keff = k + 1;

      w = L * v;
      double alpha = arma::dot(v, w);
      w -= alpha * v;
      if (k > 0) w -= betas[k - 1] * v_prev;

      double beta = arma::norm(w, 2);
      alphas.push_back(alpha);
      if (beta < 1e-14 || k == Kmax - 1) {
        break;
      }
      betas.push_back(beta);

      v_prev.swap(v);
      v = w / beta;
    }

    mat T(keff, keff, arma::fill::zeros);
    for (int k = 0; k < keff; ++k) {
      T(k, k) = alphas[k];
      if (k + 1 < keff) {
        T(k, k + 1) = betas[k];
        T(k + 1, k) = betas[k];
      }
    }

    vec eigval;
    mat eigvec;
    arma::eig_sym(eigval, eigvec, T);

    Rcpp::NumericVector f_eval_rcpp = kernel(Rcpp::wrap(eigval));
    if (f_eval_rcpp.size() != eigval.n_elem) {
      Rcpp::stop("kernel must return a vector of length equal to eigval");
    }
    vec f_eval = Rcpp::as<vec>(f_eval_rcpp);

    vec e1 = eigvec.row(0).t();
    vec fTe1 = eigvec * (f_eval % e1);
    Y.col(col) = normx * (V.head_cols(keff) * fTe1);
  }

  return Y;
}

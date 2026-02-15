#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using arma::sp_mat;
using arma::mat;
using arma::vec;

// Column-parallel sparse mat * dense mat (when OpenMP available)
// Falls back to Armadillo's built-in when OpenMP not available
static void spmm_cols(const sp_mat& A, const mat& B, mat& out) {
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  for (arma::uword j = 0; j < B.n_cols; ++j) {
    out.col(j) = A * B.col(j);
  }
#else
  out = A * B;
#endif
}

// Chebyshev polynomial filtering on a sparse Laplacian.
// coeffs: length K (K >= 1)
// X: n x m signals (columns)
// lambda_max: spectral radius estimate (>0)
// [[Rcpp::export]]
arma::mat chebyshev_filter_cpp(const arma::sp_mat& L, const arma::mat& X, const arma::vec& coeffs, double lambda_max) {
  const std::size_t K = coeffs.n_elem;
  if (K == 0) Rcpp::stop("coeffs must have length >= 1");
  if (lambda_max <= 0.0) Rcpp::stop("lambda_max must be > 0");

  const std::size_t n = L.n_rows;
  if (L.n_cols != n) Rcpp::stop("L must be square");
  if (X.n_rows != n) Rcpp::stop("X rows must match L dimension");

  const arma::uword m = X.n_cols;
  mat Y(n, m, arma::fill::zeros);
  mat T_prev = X;                  // T0
  mat T_curr(n, m);                // T1 placeholder
  mat Z(n, m);                     // scratch for L * T_k

  const double s = 2.0 / lambda_max;

  // k = 0 term
  Y += coeffs[0] * T_prev;
  if (K == 1) return Y;

  // k = 1 term: T1 = (sL - I) X = s * (L * X) - X
  spmm_cols(L, X, Z);
  T_curr = (s * Z) - X;
  Y += coeffs[1] * T_curr;
  if (K == 2) return Y;

  for (std::size_t k = 2; k < K; ++k) {
    spmm_cols(L, T_curr, Z);                    // Z = L * T_k
    Z = (2.0 * s) * Z - 2.0 * T_curr - T_prev;  // recurrence without forming Ltilde
    Y += coeffs[k] * Z;
    T_prev.swap(T_curr);
    T_curr.swap(Z);
  }

  return Y;
}

// Chebyshev filter using a cached sp_mat (external pointer) to avoid per-call conversion.
// [[Rcpp::export]]
arma::mat chebyshev_filter_xptr_cpp(SEXP L_ptr, const arma::mat& X, const arma::vec& coeffs, double lambda_max) {
    Rcpp::XPtr<arma::sp_mat> ptr(L_ptr);
    return chebyshev_filter_cpp(*ptr, X, coeffs, lambda_max);
}

// Create an external pointer to a sp_mat (for caching).
// [[Rcpp::export]]
SEXP laplacian_sp_xptr_cpp(const arma::sp_mat& L) {
    return Rcpp::XPtr<arma::sp_mat>(new arma::sp_mat(L), true);
}

// Set OpenMP threads (no-op if OpenMP unavailable)
// [[Rcpp::export]]
void set_omp_threads_cpp(int threads) {
#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#endif
}

// Chebyshev filter with optional explicit threads (double precision).
// [[Rcpp::export]]
arma::mat chebyshev_filter_omp_cpp(SEXP L_ptr,
                                   const arma::mat& X,
                                   const arma::vec& coeffs,
                                   double lambda_max,
                                   int threads = -1,
                                   int precision = 0) {
#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#endif
  Rcpp::XPtr<arma::sp_mat> ptr(L_ptr);
  const arma::sp_mat& Ld = *ptr;
  const std::size_t K = coeffs.n_elem;
  if (K == 0) Rcpp::stop("coeffs must have length >= 1");
  const std::size_t n = Ld.n_rows;
  if (X.n_rows != n) Rcpp::stop("X rows must match L dimension");

  if (precision != 0) {
    Rcpp::stop("Only double precision supported in chebyshev_filter_omp_cpp");
  }

  // Use the reference implementation (double)
  return chebyshev_filter_cpp(Ld, X, coeffs, lambda_max);
}

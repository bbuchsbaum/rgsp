#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using arma::sp_mat;
using arma::mat;
using arma::vec;

static inline void spmm_generic(const sp_mat& A, const mat& B, mat& out) {
  out = A * B;
}

template <int M>
static inline void spmm_csc_transposed_fixed(const sp_mat& A, const mat& Bt, mat& Yt) {
  const arma::uword* col_ptrs = A.col_ptrs;
  const arma::uword* row_indices = A.row_indices;
  const double* values = A.values;
  const double* x = Bt.memptr();
  double* y = Yt.memptr();

  Yt.zeros();

  for (arma::uword col = 0; col < A.n_cols; ++col) {
    const arma::uword start = col_ptrs[col];
    const arma::uword end = col_ptrs[col + 1];
    if (start == end) continue;

    const double* x_col = x + col * M;
    bool all_zero = true;
    for (int sig = 0; sig < M; ++sig) {
      all_zero = all_zero && (x_col[sig] == 0.0);
    }
    if (all_zero) {
      continue;
    }

    for (arma::uword idx = start; idx < end; ++idx) {
      const arma::uword row = row_indices[idx];
      const double val = values[idx];
      double* y_col = y + row * M;
      for (int sig = 0; sig < M; ++sig) {
        y_col[sig] += val * x_col[sig];
      }
    }
  }
}

static inline void spmm_csc_transposed_generic(const sp_mat& A, const mat& Bt, mat& Yt) {
  const arma::uword m = Bt.n_rows;
  const arma::uword* col_ptrs = A.col_ptrs;
  const arma::uword* row_indices = A.row_indices;
  const double* values = A.values;
  const double* x = Bt.memptr();
  double* y = Yt.memptr();

  Yt.zeros();

  for (arma::uword col = 0; col < A.n_cols; ++col) {
    const arma::uword start = col_ptrs[col];
    const arma::uword end = col_ptrs[col + 1];
    if (start == end) continue;

    const double* x_col = x + col * m;
    bool all_zero = true;
    for (arma::uword sig = 0; sig < m; ++sig) {
      all_zero = all_zero && (x_col[sig] == 0.0);
    }
    if (all_zero) continue;

    for (arma::uword idx = start; idx < end; ++idx) {
      const arma::uword row = row_indices[idx];
      const double val = values[idx];
      double* y_col = y + row * m;
      for (arma::uword sig = 0; sig < m; ++sig) {
        y_col[sig] += val * x_col[sig];
      }
    }
  }
}

static inline void spmm_csc_transposed(const sp_mat& A, const mat& Bt, mat& Yt) {
  switch (Bt.n_rows) {
    case 1: spmm_csc_transposed_fixed<1>(A, Bt, Yt); return;
    case 2: spmm_csc_transposed_fixed<2>(A, Bt, Yt); return;
    case 3: spmm_csc_transposed_fixed<3>(A, Bt, Yt); return;
    case 4: spmm_csc_transposed_fixed<4>(A, Bt, Yt); return;
    case 5: spmm_csc_transposed_fixed<5>(A, Bt, Yt); return;
    case 6: spmm_csc_transposed_fixed<6>(A, Bt, Yt); return;
    case 7: spmm_csc_transposed_fixed<7>(A, Bt, Yt); return;
    case 8: spmm_csc_transposed_fixed<8>(A, Bt, Yt); return;
    default: spmm_csc_transposed_generic(A, Bt, Yt); return;
  }
}

static arma::mat chebyshev_filter_smallm_cpp(const arma::sp_mat& L,
                                             const arma::mat& X,
                                             const arma::vec& coeffs,
                                             double lambda_max) {
  const arma::uword n = L.n_rows;
  const arma::uword m = X.n_cols;
  const std::size_t K = coeffs.n_elem;
  const double s = 2.0 / lambda_max;
  const double alpha = 2.0 * s;

  mat T_prev_t = X.t();
  mat T_curr_t(m, n, arma::fill::zeros);
  mat Z_t(m, n, arma::fill::zeros);
  mat Y_t = coeffs[0] * T_prev_t;
  const arma::uword len = m * n;

  if (K == 1) return Y_t.t();

  spmm_csc_transposed(L, T_prev_t, Z_t);
  {
    const double c1 = coeffs[1];
    const double* prev = T_prev_t.memptr();
    const double* z = Z_t.memptr();
    double* curr = T_curr_t.memptr();
    double* y = Y_t.memptr();
    for (arma::uword idx = 0; idx < len; ++idx) {
      curr[idx] = s * z[idx] - prev[idx];
      y[idx] += c1 * curr[idx];
    }
  }
  if (K == 2) return Y_t.t();

  for (std::size_t k = 2; k < K; ++k) {
    spmm_csc_transposed(L, T_curr_t, Z_t);
    const double ck = coeffs[k];
    double* prev = T_prev_t.memptr();
    double* curr = T_curr_t.memptr();
    double* z = Z_t.memptr();
    double* y = Y_t.memptr();
    for (arma::uword idx = 0; idx < len; ++idx) {
      const double next = alpha * z[idx] - 2.0 * curr[idx] - prev[idx];
      y[idx] += ck * next;
      prev[idx] = curr[idx];
      curr[idx] = next;
    }
  }

  return Y_t.t();
}

static inline void apply_ring_laplacian_transposed(const mat& Xt,
                                                   mat& Yt,
                                                   bool periodic) {
  const arma::uword m = Xt.n_rows;
  const arma::uword n = Xt.n_cols;
  Yt.zeros();
  if (n == 0) return;

  if (periodic) {
    for (arma::uword j = 0; j < n; ++j) {
      const arma::uword left = (j == 0) ? (n - 1) : (j - 1);
      const arma::uword right = (j + 1 == n) ? 0 : (j + 1);
      const double* xj = Xt.colptr(j);
      const double* xl = Xt.colptr(left);
      const double* xr = Xt.colptr(right);
      double* yj = Yt.colptr(j);
      for (arma::uword sig = 0; sig < m; ++sig) {
        yj[sig] = 2.0 * xj[sig] - xl[sig] - xr[sig];
      }
    }
  } else {
    if (n == 1) {
      Yt.zeros();
      return;
    }
    {
      const double* x0 = Xt.colptr(0);
      const double* x1 = Xt.colptr(1);
      double* y0 = Yt.colptr(0);
      for (arma::uword sig = 0; sig < m; ++sig) {
        y0[sig] = x0[sig] - x1[sig];
      }
    }
    for (arma::uword j = 1; j + 1 < n; ++j) {
      const double* xj = Xt.colptr(j);
      const double* xl = Xt.colptr(j - 1);
      const double* xr = Xt.colptr(j + 1);
      double* yj = Yt.colptr(j);
      for (arma::uword sig = 0; sig < m; ++sig) {
        yj[sig] = 2.0 * xj[sig] - xl[sig] - xr[sig];
      }
    }
    {
      const double* xn = Xt.colptr(n - 1);
      const double* xp = Xt.colptr(n - 2);
      double* yn = Yt.colptr(n - 1);
      for (arma::uword sig = 0; sig < m; ++sig) {
        yn[sig] = xn[sig] - xp[sig];
      }
    }
  }
}

// [[Rcpp::export]]
arma::mat chebyshev_filter_ring_cpp(const arma::mat& X,
                                    const arma::vec& coeffs,
                                    double lambda_max,
                                    bool periodic = true) {
  const std::size_t K = coeffs.n_elem;
  if (K == 0) Rcpp::stop("coeffs must have length >= 1");
  if (lambda_max <= 0.0) Rcpp::stop("lambda_max must be > 0");

  const arma::uword n = X.n_rows;
  const arma::uword m = X.n_cols;
  const double s = 2.0 / lambda_max;

  mat T_prev_t = X.t();
  mat T_curr_t(m, n, arma::fill::zeros);
  mat Z_t(m, n, arma::fill::zeros);
  mat Y_t = coeffs[0] * T_prev_t;

  if (K == 1) return Y_t.t();

  apply_ring_laplacian_transposed(T_prev_t, Z_t, periodic);
  T_curr_t = (s * Z_t) - T_prev_t;
  Y_t += coeffs[1] * T_curr_t;
  if (K == 2) return Y_t.t();

  for (std::size_t k = 2; k < K; ++k) {
    apply_ring_laplacian_transposed(T_curr_t, Z_t, periodic);
    Z_t = (2.0 * s) * Z_t - 2.0 * T_curr_t - T_prev_t;
    Y_t += coeffs[k] * Z_t;
    T_prev_t.swap(T_curr_t);
    T_curr_t.swap(Z_t);
  }

  return Y_t.t();
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
  if (m > 0 && m <= 8) {
    return chebyshev_filter_smallm_cpp(L, X, coeffs, lambda_max);
  }
  mat Y(n, m, arma::fill::zeros);
  mat T_prev = X;                  // T0
  mat T_curr(n, m);                // T1 placeholder
  mat Z(n, m);                     // scratch for L * T_k

  const double s = 2.0 / lambda_max;

  // k = 0 term
  Y += coeffs[0] * T_prev;
  if (K == 1) return Y;

  // k = 1 term: T1 = (sL - I) X = s * (L * X) - X
  spmm_generic(L, X, Z);
  T_curr = (s * Z) - X;
  Y += coeffs[1] * T_curr;
  if (K == 2) return Y;

  for (std::size_t k = 2; k < K; ++k) {
    spmm_generic(L, T_curr, Z);                    // Z = L * T_k
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

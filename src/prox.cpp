#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::vec;
using arma::mat;
using arma::sp_mat;

// Soft-threshold (L1 proximal) on a vector or each column of a matrix.
// [[Rcpp::export]]
arma::mat prox_l1_cpp(const arma::mat& X, double tau) {
  arma::mat Y = X;
  Y.transform([tau](double v) {
    double s = std::abs(v) - tau;
    return (s > 0.0) ? ((v >= 0.0 ? 1.0 : -1.0) * s) : 0.0;
  });
  return Y;
}

// Isotropic TV prox on 1D chain: prox_{lambda * TV}(x)
// Simple Chambolle projection (1997) for 1D
// [[Rcpp::export]]
arma::vec prox_tv1d_cpp(const arma::vec& y, double lambda) {
  if (lambda <= 0.0) return y;
  const arma::uword n = y.n_elem;
  if (n < 2) return y;

  arma::vec p(n - 1, arma::fill::zeros);
  arma::vec divp(n, arma::fill::zeros);
  arma::vec x(n);
  const double inv_lambda = 1.0 / lambda;
  const double tol = 1e-6;

  for (int iter = 0; iter < 100; ++iter) {
    // divergence of p
    divp[0] = -p[0];
    for (arma::uword i = 1; i < n - 1; ++i) {
      divp[i] = p[i - 1] - p[i];
    }
    divp[n - 1] = p[n - 2];

    // u = y - lambda * divp (stored in x to reuse buffer)
    x = y - lambda * divp;

    double max_delta = 0.0;
    for (arma::uword i = 0; i < n - 1; ++i) {
      double grad = x[i + 1] - x[i];
      double v = p[i] + inv_lambda * grad;
      double denom = 1.0 + inv_lambda * std::abs(grad);
      double new_p = v / denom;
      double delta = std::abs(new_p - p[i]);
      if (delta > max_delta) max_delta = delta;
      p[i] = new_p;
    }
    if (max_delta < tol) break;
  }

  // Final projection x = y - lambda * divp (divp already holds divergence of final p)
  divp[0] = -p[0];
  for (arma::uword i = 1; i < n - 1; ++i) {
    divp[i] = p[i - 1] - p[i];
  }
  divp[n - 1] = p[n - 2];
  x = y - lambda * divp;
  return x;
}

// FISTA for isotropic graph-TV denoising:
// minimize 0.5*||x - y||^2 + lambda * ||G x||_{2,1}
// L is Laplacian; we use its incidence via gradients as L times x (approx).
// Here we approximate TV with sqrt(sum_i ( (Lx)_i^2 )) which is quadratic;
// for speed we use proximal gradient with soft-threshold on Lx.
// [[Rcpp::export]]
arma::vec tv_denoise_fista_cpp(const arma::sp_mat& L, const arma::vec& y, double lambda,
                               int maxit = 200, double tol = 1e-6, double step = 1e-1) {
  arma::vec x = y;
  arma::vec z = x;
  arma::vec x_old(y.n_elem);
  arma::vec grad(y.n_elem);
  arma::vec Lz(y.n_elem);
  double t = 1.0;

  for (int k = 0; k < maxit; ++k) {
    x_old = x;

    Lz = L * z;
    grad = z;
    grad -= y;
    grad += lambda * Lz;

    x = z - step * grad;

    double t_new = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t * t));
    z = x + ((t - 1.0) / t_new) * (x - x_old);

    if (arma::norm(x - x_old, 2) < tol * std::max(1.0, arma::norm(x_old, 2))) {
      break;
    }
    t = t_new;
  }
  return x;
}

// Isotropic graph TV proximal using primal-dual (Chambolle-Pock)
// Min_x 0.5||x - y||^2 + lambda * sum_e |w_e * (x_i - x_j)|
// edges defined by rows, cols, weights w
// [[Rcpp::export]]
arma::vec prox_graph_tv_cpp(const arma::vec& y,
                            const arma::uvec& rows,
                            const arma::uvec& cols,
                            const arma::vec& w,
                            double lambda,
                            int maxit = 300,
                            double tau = 0.5,
                            double sigma = 0.5,
                            double tol = 1e-6) {
  if (lambda <= 0.0) return y;
  const arma::uword n = y.n_elem;
  const arma::uword m = rows.n_elem;

  arma::vec x = y;
  arma::vec x_bar = x;
  arma::vec x_old(n);
  arma::vec p(m, arma::fill::zeros);
  arma::vec divp(n, arma::fill::zeros);

  const arma::uword* row_ptr = rows.memptr();
  const arma::uword* col_ptr = cols.memptr();
  const double* w_ptr = w.memptr();
  double* p_ptr = p.memptr();
  double* x_bar_ptr = x_bar.memptr();

  for (int k = 0; k < maxit; ++k) {
    // dual update: p = proj_{|.|<=lambda}(p + sigma * grad x_bar)
    for (arma::uword e = 0; e < m; ++e) {
      double grad = w_ptr[e] * (x_bar_ptr[row_ptr[e]] - x_bar_ptr[col_ptr[e]]);
      double v = p_ptr[e] + sigma * grad;
      double s = std::abs(v);
      double scale = std::max(1.0, s / lambda);
      p_ptr[e] = v / scale;
    }

    // divergence of p
    divp.zeros();
    for (arma::uword e = 0; e < m; ++e) {
      double val = w_ptr[e] * p_ptr[e];
      divp[row_ptr[e]] += val;
      divp[col_ptr[e]] -= val;
    }

    // primal update: x = (y + x - tau * divp) / (1 + tau)
    x_old = x;
    x = (y + x - tau * divp) / (1.0 + tau);
    x_bar = 2.0 * x - x_old;
    x_bar_ptr = x_bar.memptr(); // refresh pointer after reallocation guard

    if (arma::norm(x - x_old, 2) < tol * std::max(1.0, arma::norm(x_old, 2))) {
      break;
    }
  }
  return x;
}

// Projection onto simplex {x >=0, sum x = 1}
// [[Rcpp::export]]
arma::vec proj_simplex_cpp(const arma::vec& y) {
    arma::vec x = y;
    arma::uword n = x.n_elem;
    arma::vec u = arma::sort(x, "descend");
    arma::vec css = arma::cumsum(u);
    arma::vec rho_vec = u - (css - 1) / arma::conv_to<arma::vec>::from(arma::regspace<arma::uvec>(1, n));
    arma::uword rho = 0;
    for (arma::uword i = 0; i < n; ++i) {
        if (rho_vec[i] > 0) rho = i;
    }
    double theta = (css[rho] - 1) / (rho + 1);
    x = arma::clamp(x - theta, 0.0, arma::datum::inf);
    return x;
}

// Projection onto l2-ball of radius r
// [[Rcpp::export]]
arma::vec proj_l2_ball_cpp(const arma::vec& y, double r) {
    double nrm = arma::norm(y, 2);
    if (nrm <= r) return y;
    return (r / nrm) * y;
}
// Box projection (component-wise clamp)
// [[Rcpp::export]]
arma::vec prox_box_cpp(const arma::vec& y, double lower, double upper) {
    arma::vec x = y;
    x.transform([&](double v) {
        if (v < lower) return lower;
        if (v > upper) return upper;
        return v;
    });
    return x;
}
// Solve (I + tau * L) X = B using conjugate gradient (symmetric positive definite).
// Returns solution matrix with same dims as B.
// [[Rcpp::export]]
arma::mat tikhonov_solve_cpp(const arma::sp_mat& L, const arma::mat& B, double tau,
                             double tol = 1e-6, int maxit = 1000) {
  arma::uword n = L.n_rows;
  mat X(B.n_rows, B.n_cols, arma::fill::zeros);
  sp_mat A = sp_mat(arma::speye(n, n)) + tau * L;

  arma::vec x(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros);
  arma::vec p_vec(n, arma::fill::zeros);
  arma::vec Ap(n, arma::fill::zeros);

  for (arma::uword c = 0; c < B.n_cols; ++c) {
    const arma::vec b = B.col(c);
    x.zeros();
    r = b;            // A * 0 = 0
    p_vec = r;
    double rsold = arma::dot(r, r);

    for (int it = 0; it < maxit; ++it) {
      Ap = A * p_vec;
      double denom = arma::dot(p_vec, Ap);
      if (denom == 0.0) break;
      double alpha = rsold / denom;
      x += alpha * p_vec;
      r -= alpha * Ap;
      double rsnew = arma::dot(r, r);
      if (std::sqrt(rsnew) < tol) break;
      p_vec = r + (rsnew / rsold) * p_vec;
      rsold = rsnew;
    }
    X.col(c) = x;
  }
  return X;
}

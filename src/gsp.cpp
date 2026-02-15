#include <RcppArmadillo.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::sp_mat;
using arma::vec;

// Compute node degrees (sum of rows) for a sparse adjacency.
// [[Rcpp::export]]
arma::vec degree_cpp(const arma::sp_mat& A) {
    arma::vec d(A.n_rows, arma::fill::zeros);
    for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
        d[it.row()] += (*it);
    }
    return d;
}

// Combinatorial or normalized Laplacian.
// [[Rcpp::export]]
arma::sp_mat laplacian_cpp(const arma::sp_mat& A, bool normalized = false, bool directed = false) {
    arma::vec out_deg = degree_cpp(A);

    if (!normalized) {
        arma::sp_mat L = -A;
        L.diag() += out_deg;
        return L;
    }

    if (!directed) {
        // Normalized: L = I - D^{-1/2} A D^{-1/2} with elementwise scaling (O(nnz)).
        arma::vec dinv_sqrt = 1.0 / arma::sqrt(out_deg);
        dinv_sqrt.replace(arma::datum::inf, 0.0);
        dinv_sqrt.replace(arma::datum::nan, 0.0);

        arma::sp_mat L = -A; // start from -A to reuse sparsity
        for (arma::sp_mat::iterator it = L.begin(); it != L.end(); ++it) {
            arma::uword i = it.row();
            arma::uword j = it.col();
            (*it) *= dinv_sqrt[i] * dinv_sqrt[j];
        }
        L.diag() += 1.0;
        return L;
    } else {
        // Random-walk normalized for directed: L = I - D_out^{-1} A (scale rows only).
        arma::vec dinv = 1.0 / out_deg;
        dinv.replace(arma::datum::inf, 0.0);
        dinv.replace(arma::datum::nan, 0.0);

        arma::sp_mat L = -A;
        for (arma::sp_mat::iterator it = L.begin(); it != L.end(); ++it) {
            arma::uword i = it.row();
            (*it) *= dinv[i];
        }
        L.diag() += 1.0;
        return L;
    }
}

// Ring/path adjacency (undirected, unweighted).
// [[Rcpp::export]]
arma::sp_mat adjacency_ring_cpp(int n, bool periodic = true) {
    if (n < 2) {
        Rcpp::stop("n must be >= 2");
    }
    int edges = periodic ? n : n - 1;
    int entries = 2 * edges; // symmetric

    arma::umat locations(2, entries);
    arma::vec values(entries, arma::fill::ones);

    int idx = 0;
    for (int i = 0; i < n - 1; ++i) {
        locations(0, idx) = i;
        locations(1, idx) = i + 1;
        ++idx;
        locations(0, idx) = i + 1;
        locations(1, idx) = i;
        ++idx;
    }
    if (periodic) {
        locations(0, idx) = 0;
        locations(1, idx) = n - 1;
        ++idx;
        locations(0, idx) = n - 1;
        locations(1, idx) = 0;
    }

    arma::sp_mat A(locations, values, n, n);
    return A;
}

// 2D grid / torus adjacency (undirected, unweighted).
// [[Rcpp::export]]
arma::sp_mat adjacency_grid2d_cpp(int nrow, int ncol, bool periodic = false) {
    if (nrow < 1 || ncol < 1) {
        Rcpp::stop("nrow and ncol must be >= 1");
    }
    const int n = nrow * ncol;
    std::vector<arma::uword> rows;
    std::vector<arma::uword> cols;
    rows.reserve(4 * n);
    cols.reserve(4 * n);

    auto id = [ncol](int r, int c) { return static_cast<arma::uword>(r * ncol + c); };
    auto add = [&](arma::uword u, arma::uword v) {
        rows.push_back(u);
        cols.push_back(v);
    };

    for (int r = 0; r < nrow; ++r) {
        for (int c = 0; c < ncol; ++c) {
            arma::uword u = id(r, c);
            int c2 = c + 1;
            if (c2 < ncol || periodic) {
                arma::uword v = id(r, (c2 >= ncol) ? 0 : c2);
                add(u, v);
                add(v, u);
            }
            int r2 = r + 1;
            if (r2 < nrow || periodic) {
                arma::uword v = id((r2 >= nrow) ? 0 : r2, c);
                add(u, v);
                add(v, u);
            }
        }
    }
    if (rows.empty()) {
        return sp_mat(n, n);
    }
    arma::umat loc(2, rows.size());
    for (size_t k = 0; k < rows.size(); ++k) {
        loc(0, k) = rows[k];
        loc(1, k) = cols[k];
    }
    arma::vec val(rows.size(), arma::fill::ones);
    return sp_mat(loc, val, n, n);
}

// 3D grid / torus adjacency (undirected, unweighted).
// [[Rcpp::export]]
arma::sp_mat adjacency_grid3d_cpp(int nx, int ny, int nz, bool periodic = false) {
    if (nx < 1 || ny < 1 || nz < 1) {
        Rcpp::stop("nx, ny, nz must be >= 1");
    }
    const int n = nx * ny * nz;
    std::vector<arma::uword> rows;
    std::vector<arma::uword> cols;
    rows.reserve(6 * n);
    cols.reserve(6 * n);

    auto id = [ny, nz](int x, int y, int z) {
        return static_cast<arma::uword>((x * ny + y) * nz + z);
    };
    auto add = [&](arma::uword u, arma::uword v) {
        rows.push_back(u);
        cols.push_back(v);
    };

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                arma::uword u = id(x, y, z);
                int x2 = x + 1;
                if (x2 < nx || periodic) {
                    arma::uword v = id((x2 >= nx) ? 0 : x2, y, z);
                    add(u, v);
                    add(v, u);
                }
                int y2 = y + 1;
                if (y2 < ny || periodic) {
                    arma::uword v = id(x, (y2 >= ny) ? 0 : y2, z);
                    add(u, v);
                    add(v, u);
                }
                int z2 = z + 1;
                if (z2 < nz || periodic) {
                    arma::uword v = id(x, y, (z2 >= nz) ? 0 : z2);
                    add(u, v);
                    add(v, u);
                }
            }
        }
    }
    if (rows.empty()) {
        return sp_mat(n, n);
    }
    arma::umat loc(2, rows.size());
    for (size_t k = 0; k < rows.size(); ++k) {
        loc(0, k) = rows[k];
        loc(1, k) = cols[k];
    }
    arma::vec val(rows.size(), arma::fill::ones);
    arma::sp_mat A(loc, val, n, n);
    return A;
}

# rgsp: Graph Signal Processing in R

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbuchsbaum/rgsp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/rgsp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**rgsp** is an R package for Graph Signal Processing (GSP), providing tools to analyze and process signals defined on graphs. It is a port of the Python [PyGSP](https://pygsp.readthedocs.io/) library with C++ acceleration via RcppArmadillo.

## Features

- **Graph Construction**: Ring, path, grid, sensor networks, random graphs (Erdos-Renyi, Barabasi-Albert, SBM)
- **Spectral Analysis**: Graph Fourier Transform (GFT/IGFT), eigendecomposition with caching
- **Filtering**: Chebyshev polynomial approximation, Lanczos, various spectral kernels
- **Wavelets**: Spectral Graph Wavelet Transform (SGWT), tight frames, multi-scale analysis
- **Optimization**: Tikhonov regularization, total variation denoising/inpainting
- **Learning**: Semi-supervised regression and classification
- **Reduction**: Kron reduction, graph multiresolution, pyramid transforms
- **Visualization**: ggplot2-based graph and signal plotting

## Installation

```r
# Install from GitHub (development version)
# install.packages("devtools")
devtools::install_github("bbuchsbaum/rgsp")
```

### Optional: enable PyGSP parity tests
- Create a local Python env with PyGSP 0.6.1 and friends:
  ```bash
  python3.12 -m venv .venv-pygsp
  . .venv-pygsp/bin/activate
  pip install numpy scipy networkx pygsp==0.6.1
  ```
- Point reticulate at it when running tests:
  ```bash
  export RETICULATE_PYTHON=$PWD/.venv-pygsp/bin/python
  ```

## Quick Start

```r
library(rgsp)

# Create a sensor network graph
g <- graph_sensor(100, seed = 42)

# Generate a smooth signal
eig <- graph_eigenpairs(g)
signal <- eig$vectors[, 2]  # Second eigenvector

# Add noise
noisy_signal <- signal + rnorm(100, sd = 0.2)

# Denoise with Tikhonov smoothing
smoothed <- tikhonov_smooth(g, noisy_signal, tau = 2)

# Or use spectral filtering
filtered <- filter_signal(g, noisy_signal, kernel_heat(1), K = 30)

# Wavelet analysis
W <- sgwt(g, signal, scales = 4)
print(W)

# Perfect reconstruction
reconstructed <- isgwt(W)
max(abs(signal - reconstructed))  # ~0
```

## Visualization

```r
library(ggplot2)

# Set up graph coordinates
g <- graph_set_coords(g, "spring", seed = 123)

# Plot graph with signal
plot_signal(g, signal)

# Plot filter responses
kernels <- list(Heat = kernel_heat(1), Mexican = kernel_mexican_hat(1))
plot_filter_response(kernels, lmax = 4)

# Wavelet coefficient heatmap
plot_wavelet_coeffs(W)
```

## PyGSP Transition Guide

| PyGSP | rgsp |
|-------|------|
| `graphs.Ring(N)` | `graph_ring(n)` |
| `graphs.Grid2d(N, M)` | `graph_grid2d(n, m)` |
| `graphs.Sensor(N)` | `graph_sensor(n)` |
| `G.compute_fourier_basis()` | `graph_eigenpairs(g)` |
| `G.gft(x)` | `gft(g, x)` |
| `G.igft(x_hat)` | `igft(g, x_hat)` |
| `filters.Heat(G, tau)` | `kernel_heat(tau)` |
| `f.filter(x)` | `filter_signal(g, x, kernel, K)` |
| `filters.MexicanHat(G, Nf)` | `mexican_hat_wavelet_bank(g, n_scales)` |
| `reduction.kron_reduction(G, keep)` | `kron_reduction(g, keep)` |
| `learning.regression_tikhonov(G, y, tau)` | `regression_tikhonov(g, y, tau)` |

## Documentation

- `vignette("getting-started")` - Introduction and basic usage
- `vignette("filters")` - Filters, wavelets, and spectral analysis
- `vignette("optimization")` - Optimization and learning algorithms

## Performance

rgsp uses RcppArmadillo for efficient sparse matrix operations. Key optimizations:

- **Chebyshev filtering**: O(K * edges) complexity, no eigendecomposition needed
- **Caching**: Eigenpairs and Chebyshev coefficients cached automatically
- **Sparse matrices**: All operations use sparse representations
- **OpenMP**: Optional parallel acceleration (when available)

Typical timings on a modern laptop (n=512 nodes):
- GFT: ~5ms
- Chebyshev filter (K=30): ~6ms
- SGWT (4 scales): ~15ms

## References

- Shuman, D. I., et al. (2013). "The emerging field of signal processing on graphs." *IEEE Signal Processing Magazine*.
- Hammond, D. K., et al. (2011). "Wavelets on graphs via spectral graph theory." *Applied and Computational Harmonic Analysis*.
- Defferrard, M., et al. (2017). "PyGSP: Graph Signal Processing in Python." https://github.com/epfl-lts2/pygsp

## License
MIT

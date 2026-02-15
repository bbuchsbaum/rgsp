# rgsp Benchmarks

This directory contains benchmark scripts for measuring rgsp performance and comparing against PyGSP.

## Benchmark Scripts

| Script | Description |
|--------|-------------|
| `bench_all.R` | Master script that runs all benchmarks |
| `bench_gft.R` | GFT/IGFT performance and comparison with PyGSP |
| `bench_chebyshev.R` | Chebyshev polynomial filtering benchmarks |
| `bench_wavelets.R` | SGWT wavelet transform benchmarks |
| `bench_filters.R` | General filter bank benchmarks |
| `bench_threads.R` | OpenMP/BLAS thread scaling benchmarks |
| `bench_tv.R` | Total variation optimization benchmarks |
| `bench_pygsp_compare.R` | Side-by-side timing vs PyGSP (GFT/IGFT, heat filter) |
| `bench_sparsify_spectral.R` | Spectral sparsification + downstream speedup benchmark |

## Running Benchmarks

### Quick Smoke Test (CI)

```r
Rscript inst/bench/bench_all.R --quick
```

Runs minimal benchmarks to verify performance hasn't regressed. Completes in ~30 seconds.

### Full Benchmark Suite

```r
Rscript inst/bench/bench_all.R
```

Runs comprehensive benchmarks across multiple graph sizes, filter orders, and signal counts. Requires `bench` package and optionally `reticulate` with PyGSP for comparisons.

### Individual Benchmarks

```r
# GFT/IGFT benchmarks
Rscript inst/bench/bench_gft.R

# Chebyshev filtering benchmarks
Rscript inst/bench/bench_chebyshev.R

# Wavelet transform benchmarks
Rscript inst/bench/bench_wavelets.R

# Thread scaling benchmarks
Rscript inst/bench/bench_threads.R
```

## Comparing with PyGSP

To enable PyGSP comparisons, install PyGSP in a Python environment accessible to reticulate:

```bash
pip install pygsp numpy
```

Example one-off comparison (writes CSV when `BENCH_OUT` set):

```bash
RETICULATE_PYTHON=.venv-pygsp/bin/python \
BENCH_OUT=bench_pygsp.csv \
Rscript inst/bench/bench_pygsp_compare.R
```

The benchmark scripts automatically detect PyGSP availability and include comparisons when possible.

## Typical Results

### GFT/IGFT (Ring graph, n=512)

| Operation | rgsp (ms) | PyGSP (ms) | Speedup |
|-----------|-----------|------------|---------|
| Eigenpairs | ~50 | ~80 | 1.6x |
| GFT | ~2 | ~3 | 1.5x |
| IGFT | ~2 | ~3 | 1.5x |

### Chebyshev Filtering (K=30, 4 signals)

| Graph Size | rgsp (ms) | PyGSP (ms) | Speedup |
|------------|-----------|------------|---------|
| n=128 | ~1.5 | ~3 | 2x |
| n=512 | ~6 | ~15 | 2.5x |
| n=2048 | ~45 | ~120 | 2.7x |

### SGWT (4 scales)

| Graph Size | rgsp Analysis (ms) | rgsp Synthesis (ms) |
|------------|-------------------|---------------------|
| n=128 | ~4 | ~5 |
| n=512 | ~15 | ~18 |
| n=1024 | ~35 | ~40 |

*Results vary based on hardware, BLAS implementation, and compiler optimizations.*

## Performance Notes

### BLAS Provider

rgsp uses RcppArmadillo which leverages your system's BLAS implementation. For best performance:

- **macOS**: Uses Accelerate framework by default (excellent performance)
- **Linux**: Install OpenBLAS or Intel MKL
- **Windows**: R's reference BLAS is slower; consider Microsoft R Open with MKL

Check your BLAS with:

```r
sessionInfo()  # Look for BLAS/LAPACK info
```

### OpenMP

The C++ backend supports OpenMP for parallel operations when available. Thread count can be controlled via:

```r
# If RhpcBLASctl is installed
RhpcBLASctl::omp_set_num_threads(4)
```

### Caching

rgsp caches:
- Graph eigenpairs (for GFT/IGFT)
- Chebyshev coefficients (for filtering)

Clear caches with:

```r
graph_cache_clear()
cheby_cache_clear()
```

## CI Integration

The GitHub Actions workflow `.github/workflows/benchmarks.yaml` runs quick benchmarks on releases and can be manually triggered for full benchmarks.

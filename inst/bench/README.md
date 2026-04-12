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
| `bench_pygsp_compare.R` | Side-by-side timing vs PyGSP for ring/path/grid structured cases |
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
BENCH_OUT=bench_pygsp.csv \
Rscript inst/bench/bench_pygsp_compare.R
```

`bench_pygsp_compare.R` loads the local working tree with `pkgload::load_all()` and
automatically uses the vendored `./pygsp` when present, so it measures the code you
are editing rather than an older installed package.

## Interpreting Results

Benchmark results are highly hardware-dependent and are sensitive to:

- the active BLAS/LAPACK implementation used by R,
- the Python and SciPy wheels backing PyGSP,
- whether `reticulate` is pointed at the vendored `./pygsp` tree,
- graph size, signal count, and approximation order.

For that reason, this directory no longer publishes static "typical" speed tables.
Use `bench_pygsp_compare.R`, `bench_gft.R`, and `bench_chebyshev.R` to generate
measured results on the machine and toolchain you actually care about.

The PyGSP comparison script now includes structured heat-filter cases for:

- `ring`, where `rgsp` can use an FFT-based exact path,
- `path`, where `rgsp` can use a DCT/FFT-based exact path,
- non-periodic `grid2d`, where `rgsp` can use a separable exact path,
- generic `sensor`, `random_geometric`, and `sbm` graphs, all built from the exact same sparse adjacency in both libraries.

The benchmark output now includes `rgsp_mad_ms` and `pygsp_mad_ms` columns so
you can distinguish a real speed gap from timer noise.

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

# rgsp 0.1.0

Initial release of rgsp, a Graph Signal Processing toolkit for R.

## Features

### Graph Construction
- Basic graphs: `graph_ring()`, `graph_path()`, `graph_grid2d()`, `graph_grid3d()`
- Random graphs: `graph_erdos_renyi()`, `graph_barabasi_albert()`, `graph_random_regular()`, `graph_random_geometric()`
- Structured graphs: `graph_sensor()`, `graph_sbm()`, `graph_complete()`
- Graph properties: `graph_laplacian()`, `graph_degree()`, `graph_incidence()`
- Utilities: `graph_rescale()`, `graph_sparsify()`, `lambda_max()`

### Spectral Analysis
- Graph Fourier Transform: `gft()`, `igft()`
- Eigendecomposition with caching: `graph_eigenpairs()`
- Cache management: `graph_cache_clear()`, `graph_cache_invalidate_spectra()`

### Filtering
- Chebyshev polynomial filtering: `filter_signal()`, `chebyshev_filter()`
- Lanczos filtering: `filter_signal_lanczos()`
- Filter banks: `filter_bank_apply()`
- Spectral kernels: `kernel_heat()`, `kernel_exponential()`, `kernel_mexican_hat()`, `kernel_meyer()`, `kernel_rectangle()`, `kernel_half_cosine()`, `kernel_gabor()`, `kernel_wave()`, `kernel_tight_band()`
- Heat propagation: `heat_propagate()`

### Wavelets
- SGWT: `sgwt()`, `isgwt()`
- Wavelet banks: `mexican_hat_wavelet_bank()`, `meyer_wavelet_bank()`, `wavelet_heat_bank()`, `wavelet_heat_transform()`
- Tight frames: `tight_frame_regular()`, `tight_frame_held()`, `tight_frame_simoncelli()`, `tight_frame_papadakis()`, `tight_frame_bank()`
- Gabor/modulation: `gabor_filter_bank()`, `modulation_transform()`
- Analysis tools: `sgwt_coeffs_at_scale()`, `sgwt_coeffs_tidy()`, `band_energy_df()`, `sgwt_atom()`, `sgwt_node_contribution()`, `sgwt_reconstruct_at()`, `wavelet_stat_map()`
- Frame analysis: `compute_frame()`, `compute_spectrogram()`

### Optimization
- Tikhonov smoothing: `tikhonov_smooth()`
- Total variation: `tv_denoise()`, `tv_inpaint()`
- Interpolation: `interpolate_laplacian()`
- Proximal operators: `prox_l1()`, `prox_box()`, `prox_tv1d()`, `proj_simplex()`, `proj_l2_ball()`

### Learning
- Semi-supervised learning: `regression_tikhonov()`, `classification_tikhonov()`

### Reduction
- Kron reduction: `kron_reduction()`
- Multiresolution: `graph_multiresolution()`, `graph_compress()`
- Pyramid transforms: `pyramid_analysis()`, `pyramid_synthesis()`

### Differential Operators
- Gradient and divergence: `grad()`, `div()`, `compute_differential_operator()`

### Visualization (requires ggplot2)
- Layout algorithms: `graph_set_coords()` (spring, random, circle, grid, spectral, line)
- Plotting: `plot_graph()`, `plot_signal()`, `plot_filter_response()`, `plot_wavelet_coeffs()`, `plot_band_energy()`, `plot_spectrogram()`
- Style presets: `theme_gsp()`, `scale_signal()`, `scale_fill_gsp()`, `legend_position()`, `graph_title()`

## Technical Details

- C++ backend via RcppArmadillo for efficient sparse matrix operations
- Optional OpenMP acceleration
- Chebyshev coefficient caching for repeated filtering
- Eigenpair caching with automatic invalidation
- Jackson damping for improved Chebyshev convergence

## Compatibility

- R >= 4.0
- Tested on macOS, Linux, Windows

## Acknowledgments

This package is a port of PyGSP (https://github.com/epfl-lts2/pygsp) by EPFL LTS2.

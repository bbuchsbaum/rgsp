# Graph Signal Processing Visualization
#
# Provides layout algorithms and ggplot2-based plotting for graphs and signals.
# Plotting functions require ggplot2 to be installed (suggested dependency).

# ==============================================================================
# Layout Algorithms
# ==============================================================================

#' Set or Compute Graph Coordinates
#'
#' Compute node coordinates for graph visualization using various layout algorithms.
#'
#' @param g gsp_graph object
#' @param kind layout algorithm or coordinates:
#'   - "spring" (default): Fruchterman-Reingold force-directed layout
#'   - "random": random 2D coordinates
#'   - "circle" or "ring": nodes arranged in a circle
#'   - "grid": grid layout (for grid graphs)
#'   - "spectral": Laplacian eigenmap (uses 2nd and 3rd eigenvectors
#'   - numeric matrix: custom coordinates (n x 2 or n x 3)
#' @param dim integer; 2 or 3 for 2D/3D layouts (default: 2)
#' @param seed random seed for reproducible layouts
#' @param ... additional arguments passed to layout algorithms
#'
#' @return gsp_graph with updated coords field
#' @export
#'
#' @examples
#' g <- graph_erdos_renyi(30, 0.2, seed = 42)
#' g <- graph_set_coords(g, "spring", seed = 123)
#' g$coords[1:5, ]
#'
#' # Spectral layout uses graph structure
#' g <- graph_set_coords(g, "spectral")
#'
graph_set_coords <- function(g, kind = "spring", dim = 2, seed = NULL, ...) {
  if (!inherits(g, "gsp_graph")) {
    stop("'g' must be a gsp_graph object")
  }

  n <- g$n

  # Handle numeric coordinates directly

if (is.numeric(kind) && !is.null(dim(kind))) {
    coords <- as.matrix(kind)
    if (nrow(coords) != n) {
      stop("Coordinate matrix must have ", n, " rows (one per node)")
    }
    if (ncol(coords) < 2 || ncol(coords) > 3) {
      stop("Coordinates must be 2D or 3D (2 or 3 columns)")
    }
    g$coords <- coords
    return(g)
  }

  # Set seed if provided
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    set.seed(seed)
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  coords <- switch(kind,
    spring = layout_spring(g, dim = dim, ...),
    random = layout_random(n, dim = dim),
    circle = layout_circle(n),
    ring = layout_circle(n),
    grid = layout_grid(g),
    line = layout_line(n),
    spectral = layout_spectral(g, dim = dim),
    stop("Unknown layout kind: ", kind)
  )

  g$coords <- coords
  g
}


#' Random Layout
#' @param n number of nodes
#' @param dim 2 or 3
#' @return matrix of coordinates
#' @keywords internal
layout_random <- function(n, dim = 2) {
  matrix(runif(n * dim), nrow = n, ncol = dim)
}


#' Circle/Ring Layout
#' @param n number of nodes
#' @return matrix of 2D coordinates
#' @keywords internal
layout_circle <- function(n) {
  angles <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  cbind(cos(angles), sin(angles))
}


#' Line Layout
#' @param n number of nodes
#' @return matrix of 2D coordinates
#' @keywords internal
layout_line <- function(n) {
  cbind(seq_len(n), rep(0, n))
}


#' Grid Layout (for grid graphs)
#' @param g graph object
#' @return matrix of coordinates
#' @keywords internal
layout_grid <- function(g) {
  # If coords already exist (from grid construction), use them
 if (!is.null(g$coords)) {
    return(g$coords)
  }
  # Otherwise, try to infer grid dimensions
  n <- g$n
  side <- floor(sqrt(n))
  cbind(
    rep(seq_len(side), length.out = n),
    rep(seq_len(ceiling(n / side)), each = side)[seq_len(n)]
  )
}


#' Spectral Layout (Laplacian Eigenmap)
#' @param g graph object
#' @param dim 2 or 3
#' @return matrix of coordinates
#' @keywords internal
layout_spectral <- function(g, dim = 2) {
  # Use 2nd through (dim+1)th eigenvectors
  eig <- graph_eigenpairs(g, k = dim + 1)
  # Skip first eigenvector (constant for connected graphs)
  eig$vectors[, 2:(dim + 1), drop = FALSE]
}


#' Fruchterman-Reingold Spring Layout
#'
#' Force-directed graph layout using the Fruchterman-Reingold algorithm.
#'
#' @param g graph object
#' @param dim 2 or 3
#' @param iterations number of iterations (default: 50)
#' @param k optimal distance between nodes (default: auto)
#' @param temp initial temperature (default: 0.1)
#' @return matrix of coordinates
#' @keywords internal
layout_spring <- function(g, dim = 2, iterations = 50, k = NULL, temp = 0.1) {
  n <- g$n
  A <- g$adjacency

  # Initial random positions
  pos <- matrix(runif(n * dim), nrow = n, ncol = dim)

  # Optimal distance
  if (is.null(k)) {
    k <- sqrt(1.0 / n)
  }

  # Convert to list-of-lists for efficient neighbor access
  A_lil <- as(A, "TsparseMatrix")

  dt <- temp / (iterations + 1)

  for (iter in seq_len(iterations)) {
    displacement <- matrix(0, nrow = n, ncol = dim)

    for (i in seq_len(n)) {
      # Difference between this node and all others
      delta <- sweep(pos, 2, pos[i, ], "-")  # pos - pos[i, ]
      delta <- -delta  # pos[i, ] - pos

      # Distance
      distance <- sqrt(rowSums(delta^2))
      distance[distance < 0.01] <- 0.01

      # Repulsive force (all pairs)
      repulsion <- delta * (k^2 / distance^2)

      # Attractive force (neighbors only)
      neighbors <- which(A[i, ] != 0)
      attraction <- matrix(0, nrow = n, ncol = dim)
      if (length(neighbors) > 0) {
        attraction[neighbors, ] <- -delta[neighbors, , drop = FALSE] *
          (distance[neighbors] / k)
      }

      # Net displacement for node i
      displacement[i, ] <- colSums(repulsion) + colSums(attraction)
    }

    # Update positions with temperature limiting
    length_disp <- sqrt(rowSums(displacement^2))
    length_disp[length_disp < 0.01] <- 0.1
    pos <- pos + displacement * (temp / length_disp)

    # Cool down
    temp <- temp - dt
  }

  # Center and scale
  pos <- scale(pos, center = TRUE, scale = FALSE)
  max_coord <- max(abs(pos))
  if (max_coord > 0) {
    pos <- pos / max_coord
  }

  pos
}


# ==============================================================================
# Plotting Functions (ggplot2-based)
# ==============================================================================

#' Check if ggplot2 is available
#' @keywords internal
check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. ",
         "Install it with: install.packages('ggplot2')")
  }
}


#' Plot a Graph
#'
#' Visualize a graph using ggplot2, optionally showing a signal as node color.
#'
#' @param g gsp_graph object
#' @param signal optional numeric vector to show as node color
#' @param vertex_size node size (default: 3)
#' @param vertex_color node color when no signal (default: "steelblue")
#' @param edge_color edge color (default: "gray70")
#' @param edge_width edge width (default: 0.3)
#' @param edge_alpha edge transparency (default: 0.5)
#' @param show_edges logical; draw edges? (default: TRUE if < 5000 edges)
#' @param layout layout algorithm if coords not set (default: "spring")
#' @param colorbar logical; show colorbar for signal (default: TRUE)
#' @param limits color scale limits c(min, max) (default: signal range)
#' @param palette color palette for signal (default: "viridis")
#' @param title plot title
#' @param seed random seed for layout
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_sensor(50, seed = 42)
#' plot_graph(g)
#'
#' # With signal as color
#' signal <- sin(seq(0, 2*pi, length.out = 50))
#' plot_graph(g, signal = signal, palette = "RdBu")
#'
#' # Custom layout
#' g <- graph_ring(20)
#' plot_graph(g, layout = "circle")
#' }
#'
plot_graph <- function(g,
                       signal = NULL,
                       vertex_size = 3,
                       vertex_color = "steelblue",
                       edge_color = "gray70",
                       edge_width = 0.3,
                       edge_alpha = 0.5,
                       show_edges = NULL,
                       layout = "spring",
                       colorbar = TRUE,
                       limits = NULL,
                       palette = "viridis",
                       title = NULL,
                       seed = NULL) {
  check_ggplot2()

  if (!inherits(g, "gsp_graph")) {
    stop("'g' must be a gsp_graph object")
  }

  # Ensure coordinates exist
 if (is.null(g$coords)) {
    g <- graph_set_coords(g, kind = layout, seed = seed)
  }

  coords <- g$coords
  n <- g$n

  # Determine whether to show edges
  n_edges <- sum(g$adjacency != 0) / 2
  if (is.null(show_edges)) {
    show_edges <- n_edges < 5000
  }

  # Create node data frame
  node_df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    id = seq_len(n)
  )

  if (!is.null(signal)) {
    if (length(signal) != n) {
      stop("signal must have length ", n)
    }
    node_df$signal <- signal
  }

  # Create edge data frame if needed
  edge_df <- NULL
  if (show_edges) {
    A <- as(g$adjacency, "TsparseMatrix")
    # Only upper triangle for undirected
    mask <- A@i < A@j
    sources <- A@i[mask] + 1L
    targets <- A@j[mask] + 1L

    if (length(sources) > 0) {
      edge_df <- data.frame(
        x = coords[sources, 1],
        y = coords[sources, 2],
        xend = coords[targets, 1],
        yend = coords[targets, 2]
      )
    }
  }

  # Build plot
  p <- ggplot2::ggplot()

  # Add edges
 if (!is.null(edge_df) && nrow(edge_df) > 0) {
    p <- p + ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = edge_color,
      linewidth = edge_width,
      alpha = edge_alpha
    )
  }

  # Add nodes
  if (!is.null(signal)) {
    p <- p + ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$signal),
      size = vertex_size
    )

    # Color scale
    if (palette == "viridis") {
      p <- p + ggplot2::scale_color_viridis_c(limits = limits)
    } else {
      p <- p + ggplot2::scale_color_distiller(
        palette = palette,
        limits = limits,
        direction = 1
      )
    }

    if (!colorbar) {
      p <- p + ggplot2::guides(color = "none")
    }
  } else {
    p <- p + ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = vertex_color,
      size = vertex_size
    )
  }

  # Theme
  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  p
}


#' Plot a Graph Signal
#'
#' Quick visualization of a signal on a graph. Wrapper around plot_graph.
#'
#' @param g gsp_graph object
#' @param signal numeric vector of signal values
#' @param ... additional arguments passed to plot_graph
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_sensor(100, seed = 42)
#' signal <- rnorm(100)
#' plot_signal(g, signal)
#' }
#'
plot_signal <- function(g, signal, ...) {
  plot_graph(g, signal = signal, ...)
}


#' Plot Filter Response
#'
#' Visualize the spectral response of filter kernels.
#'
#' @param kernels list of kernel functions, or a single kernel
#' @param lmax maximum eigenvalue (spectral radius)
#' @param n number of points for evaluation (default: 200)
#' @param eigenvalues optional vector of graph eigenvalues to mark
#' @param show_sum logical; show sum of squared responses (default: TRUE if multiple kernels)
#' @param labels logical; add legend labels (default: TRUE if multiple kernels)
#' @param title plot title
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Single kernel
#' plot_filter_response(kernel_heat(1), lmax = 4)
#'
#' # Filter bank
#' g <- graph_ring(20)
#' bank <- mexican_hat_wavelet_bank(g, n_scales = 4)
#' plot_filter_response(bank, lmax = lambda_max(g), show_sum = TRUE)
#' }
#'
plot_filter_response <- function(kernels,
                                  lmax,
                                  n = 200,
                                  eigenvalues = NULL,
                                  show_sum = NULL,
                                  labels = NULL,
                                  title = "Filter Response") {
  check_ggplot2()

  # Ensure kernels is a list
  if (is.function(kernels)) {
    kernels <- list(kernels)
  }

  n_filters <- length(kernels)

  if (is.null(show_sum)) {
    show_sum <- n_filters > 1
  }
  if (is.null(labels)) {
    labels <- n_filters > 1
  }

  # Evaluate filters
  lambda <- seq(0, lmax, length.out = n)
  responses <- sapply(kernels, function(k) k(lambda))
  if (n_filters == 1) {
    responses <- matrix(responses, ncol = 1)
  }

  # Build data frame for plotting
  df <- data.frame(
    lambda = rep(lambda, n_filters),
    response = as.vector(responses),
    filter = rep(paste0("g", seq_len(n_filters)), each = n)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lambda, y = .data$response))

  if (labels) {
    p <- p + ggplot2::geom_line(ggplot2::aes(color = .data$filter), alpha = 0.7)
  } else {
    p <- p + ggplot2::geom_line(ggplot2::aes(group = .data$filter), alpha = 0.7,
                                 color = "steelblue")
  }

  # Add sum of squares
  if (show_sum) {
    sum_sq <- rowSums(responses^2)
    sum_df <- data.frame(lambda = lambda, sum_sq = sum_sq)
    p <- p + ggplot2::geom_line(
      data = sum_df,
      ggplot2::aes(x = .data$lambda, y = .data$sum_sq),
      color = "black",
      linewidth = 1,
      linetype = "dashed"
    )
  }

  # Add eigenvalue markers if provided
 if (!is.null(eigenvalues)) {
    eig_df <- data.frame(lambda = eigenvalues)
    p <- p + ggplot2::geom_vline(
      data = eig_df,
      ggplot2::aes(xintercept = .data$lambda),
      color = "gray80",
      alpha = 0.5
    )
  }

  p <- p +
    ggplot2::labs(
      x = expression(lambda),
      y = expression(g(lambda)),
      title = title
    ) +
    ggplot2::theme_minimal()

  if (!labels) {
    p <- p + ggplot2::guides(color = "none")
  }

  p
}


#' Plot Wavelet Coefficients Heatmap
#'
#' Visualize SGWT coefficients as a heatmap (nodes x scales).
#'
#' @param W sgwt_coeffs object from sgwt()
#' @param signal_idx which signal to plot if multiple (default: 1)
#' @param title plot title
#' @param palette color palette (default: "RdBu")
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_sensor(50, seed = 42)
#' W <- sgwt(g, rnorm(50))
#' plot_wavelet_coeffs(W)
#' }
#'
plot_wavelet_coeffs <- function(W, signal_idx = 1, title = "Wavelet Coefficients",
                                 palette = "RdBu") {
  check_ggplot2()

  if (!inherits(W, "sgwt_coeffs")) {
    stop("'W' must be an sgwt_coeffs object from sgwt()")
  }

  coeffs <- W$coeffs
  dims <- dim(coeffs)

  if (length(dims) == 3) {
    coeffs_mat <- coeffs[, , signal_idx]
  } else {
    coeffs_mat <- coeffs[, , 1]
  }

  n_nodes <- nrow(coeffs_mat)
  n_bands <- ncol(coeffs_mat)

  band_names <- if (W$include_lowpass) {
    c("LP", paste0("s=", round(W$scales, 1)))
  } else {
    paste0("s=", round(W$scales, 1))
  }

  # Create data frame
  df <- expand.grid(node = seq_len(n_nodes), band = seq_len(n_bands))
  df$coeff <- as.vector(coeffs_mat)
  df$band_name <- factor(band_names[df$band], levels = band_names)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$band_name, y = .data$node,
                                         fill = .data$coeff)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_distiller(palette = palette, direction = 1) +
    ggplot2::labs(x = "Scale", y = "Node", fill = "Coeff", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}


#' Plot Band Energy Distribution
#'
#' Bar plot showing energy distribution across wavelet scales.
#'
#' @param W sgwt_coeffs object from sgwt()
#' @param title plot title
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_sensor(50, seed = 42)
#' W <- sgwt(g, rnorm(50))
#' plot_band_energy(W)
#' }
#'
plot_band_energy <- function(W, title = "Wavelet Energy by Scale") {
  check_ggplot2()

  if (!inherits(W, "sgwt_coeffs")) {
    stop("'W' must be an sgwt_coeffs object from sgwt()")
  }

  df <- band_energy_df(W)
  df$band <- factor(df$band, levels = df$band)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$band, y = .data$pct)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", .data$pct)),
                       vjust = -0.5, size = 3) +
    ggplot2::labs(x = "Scale", y = "Energy (%)", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ylim(0, max(df$pct) * 1.15)

  p
}


#' Plot Graph Spectrogram
#'
#' Visualize the vertex-frequency representation of a signal.
#'
#' @param g gsp_graph object
#' @param signal numeric vector
#' @param scales scales for the spectrogram
#' @param title plot title
#' @param palette color palette (default: "viridis")
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_ring(30)
#' signal <- sin(seq(0, 4*pi, length.out = 30))
#' plot_spectrogram(g, signal, scales = c(0.5, 1, 2, 4, 8))
#' }
#'
plot_spectrogram <- function(g, signal, scales = NULL, title = "Graph Spectrogram",
                              palette = "viridis") {
  check_ggplot2()

  if (is.null(scales)) {
    lmax <- lambda_max(g)
    scales <- lmax / 2^(0:4)
  }

  spec <- compute_spectrogram(g, signal, scales)

  n_nodes <- nrow(spec)
  n_scales <- ncol(spec)

  df <- expand.grid(node = seq_len(n_nodes), scale_idx = seq_len(n_scales))
  df$energy <- as.vector(spec)
  df$scale <- scales[df$scale_idx]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$scale), y = .data$node,
                                         fill = .data$energy)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(x = "Scale", y = "Node", fill = "Energy", title = title) +
    ggplot2::theme_minimal()

  p
}


# ==============================================================================
# Style Presets and Helpers
# ==============================================================================

#' GSP Plotting Theme Presets
#'
#' Pre-defined ggplot2 themes for graph signal processing visualizations.
#'
#' @param preset one of "default", "minimal", "dark", "publication"
#' @param base_size base font size
#'
#' @return ggplot2 theme object
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_ring(20)
#' g <- graph_set_coords(g, "circle")
#' plot_graph(g) + theme_gsp("minimal")
#' }
#'
theme_gsp <- function(preset = "default", base_size = 11) {
  check_ggplot2()

  switch(preset,
    default = ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      ),
    minimal = ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      ),
    dark = ggplot2::theme_dark(base_size = base_size) +
      ggplot2::theme(
        panel.grid = ggplot2::element_line(colour = "grey30"),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      ),
    publication = ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(fill = NA, colour = "black"),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        legend.position = "right",
        axis.line = ggplot2::element_blank()
      ),
    stop("Unknown preset: ", preset, ". Options: default, minimal, dark, publication")
  )
}


#' Color Scale Presets for Graph Signals
#'
#' Get pre-defined color scales for signal visualization.
#'
#' @param palette one of "viridis", "spectral", "diverging", "sequential", "heat"
#' @param direction 1 for normal direction, -1 for reversed
#' @param midpoint for diverging scales, the value at the midpoint (default: 0)
#'
#' @return ggplot2 scale object
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_ring(20)
#' g <- graph_set_coords(g, "circle")
#' signal <- sin(seq(0, 2*pi, length.out = 20))
#' plot_signal(g, signal) + scale_signal("diverging")
#' }
#'
scale_signal <- function(palette = "viridis", direction = 1, midpoint = 0) {
  check_ggplot2()

  switch(palette,
    viridis = ggplot2::scale_colour_viridis_c(direction = direction),
    spectral = ggplot2::scale_colour_distiller(palette = "Spectral", direction = direction),
    diverging = ggplot2::scale_colour_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = midpoint
    ),
    sequential = ggplot2::scale_colour_distiller(palette = "YlOrRd", direction = direction),
    heat = ggplot2::scale_colour_gradientn(
      colours = c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred")
    ),
    stop("Unknown palette: ", palette, ". Options: viridis, spectral, diverging, sequential, heat")
  )
}


#' Fill Scale Presets for Heatmaps
#'
#' Get pre-defined fill scales for spectrogram and coefficient visualizations.
#'
#' @param palette one of "viridis", "spectral", "magma", "inferno", "plasma"
#' @param direction 1 for normal direction, -1 for reversed
#'
#' @return ggplot2 scale object
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_ring(12)
#' W <- sgwt(g, rnorm(12))
#' plot_wavelet_coeffs(W) + scale_fill_gsp("magma")
#' }
#'
scale_fill_gsp <- function(palette = "viridis", direction = 1) {
  check_ggplot2()

  switch(palette,
    viridis = ggplot2::scale_fill_viridis_c(direction = direction),
    spectral = ggplot2::scale_fill_distiller(palette = "Spectral", direction = direction),
    magma = ggplot2::scale_fill_viridis_c(option = "magma", direction = direction),
    inferno = ggplot2::scale_fill_viridis_c(option = "inferno", direction = direction),
    plasma = ggplot2::scale_fill_viridis_c(option = "plasma", direction = direction),
    stop("Unknown palette: ", palette, ". Options: viridis, spectral, magma, inferno, plasma")
  )
}


#' Legend Position Helper
#'
#' Convenience function to adjust legend position.
#'
#' @param position one of "right", "left", "top", "bottom", "none"
#'
#' @return ggplot2 theme element for legend position
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_ring(20)
#' g <- graph_set_coords(g, "circle")
#' plot_graph(g, signal = rnorm(20)) + legend_position("bottom")
#' }
#'
legend_position <- function(position = "right") {
  check_ggplot2()
  ggplot2::theme(legend.position = position)
}


#' Add Graph Title and Subtitle
#'
#' Convenience function to add title and subtitle to graph plots.
#'
#' @param title main title
#' @param subtitle optional subtitle
#' @param caption optional caption (appears at bottom)
#'
#' @return ggplot2 labs object
#' @export
#'
#' @examples
#' \dontrun{
#' g <- graph_sensor(50, seed = 42)
#' g <- graph_set_coords(g, "spring")
#' plot_graph(g) + graph_title("Sensor Network", "n = 50 nodes")
#' }
#'
graph_title <- function(title, subtitle = NULL, caption = NULL) {
  check_ggplot2()
  ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
}

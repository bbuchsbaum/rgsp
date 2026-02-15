# Graph interoperability (Python: NetworkX / graph-tool) via reticulate

.require_reticulate <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("This function requires the 'reticulate' package.", call. = FALSE)
  }
}

.require_py_module <- function(module) {
  .require_reticulate()
  if (!reticulate::py_available(initialize = TRUE)) {
    stop("Python is not available for reticulate.", call. = FALSE)
  }
  if (!reticulate::py_module_available(module)) {
    stop("Python module '", module, "' is not available.", call. = FALSE)
  }
}

.as_scipy_csc <- function(A) {
  if (!inherits(A, "dgCMatrix")) {
    A <- methods::as(A, "dgCMatrix")
  }
  A <- Matrix::drop0(A)

  np <- reticulate::import("numpy", delay_load = FALSE)
  sp <- reticulate::import("scipy.sparse", delay_load = FALSE)

  data <- np$array(as.numeric(A@x))
  indices <- np$array(as.integer(A@i))
  indptr <- np$array(as.integer(A@p))
  shape <- reticulate::tuple(as.integer(dim(A)))

  sp$csc_matrix(reticulate::tuple(list(data, indices, indptr)), shape = shape)
}

.dgC_from_scipy_csc <- function(S) {
  shape <- as.integer(S$shape)
  A <- methods::new("dgCMatrix")
  A@Dim <- shape
  A@x <- as.numeric(S$data)
  A@i <- as.integer(S$indices)
  A@p <- as.integer(S$indptr)
  Matrix::drop0(A)
}

#' Export a graph to NetworkX
#'
#' Converts an `rgsp` graph to a Python NetworkX graph via `reticulate`.
#' Edge weights are stored under the standard NetworkX edge attribute `"weight"`.
#'
#' @param g A `gsp_graph`.
#' @param include_coords If `TRUE` and `g$coords` is present, store coordinates
#'   as node attribute `"coords"`.
#' @return A NetworkX graph object (Python).
#' @export
graph_to_networkx <- function(g, include_coords = TRUE) {
  if (!inherits(g, "gsp_graph")) stop("'g' must be a gsp_graph", call. = FALSE)
  .require_py_module("networkx")
  .require_py_module("scipy.sparse")

  nx <- reticulate::import("networkx", delay_load = FALSE)

  S <- .as_scipy_csc(g$adjacency)
  create_using <- if (isTRUE(g$directed)) nx$DiGraph() else nx$Graph()

  if (!is.null(nx$from_scipy_sparse_array)) {
    g_nx <- nx$from_scipy_sparse_array(S, create_using = create_using)
  } else {
    g_nx <- nx$from_scipy_sparse_matrix(S, create_using = create_using)
  }

  if (isTRUE(include_coords) && !is.null(g$coords)) {
    coords <- as.matrix(g$coords)
    attr_list <- lapply(seq_len(nrow(coords)), function(i) as.numeric(coords[i, ]))
    names(attr_list) <- as.character(seq_len(nrow(coords)) - 1L)
    nx$set_node_attributes(g_nx, attr_list, name = "coords")
  }

  g_nx
}

#' Import a graph from NetworkX
#'
#' Converts a Python NetworkX graph to an `rgsp` `gsp_graph` via `reticulate`.
#'
#' @param graph A NetworkX graph object (Python).
#' @param weight Edge attribute name to use as weights (default `"weight"`).
#' @param normalized Logical; whether the returned graph uses a normalized Laplacian.
#' @param include_coords If `TRUE`, attempt to import node attribute `"coords"`
#'   into `g$coords` (if present).
#' @return A `gsp_graph`.
#' @export
graph_from_networkx <- function(graph,
                                weight = "weight",
                                normalized = FALSE,
                                include_coords = TRUE) {
  .require_py_module("networkx")
  .require_py_module("scipy.sparse")

  nx <- reticulate::import("networkx", delay_load = FALSE)

  to_sp <- nx$to_scipy_sparse_array %||% nx$to_scipy_sparse_matrix
  if (is.null(to_sp)) stop("NetworkX does not expose to_scipy_sparse_* helpers.", call. = FALSE)

  S <- to_sp(graph, weight = weight, format = "csc")
  A <- .dgC_from_scipy_csc(S)

  directed <- isTRUE(graph$is_directed())

  coords <- NULL
  if (isTRUE(include_coords)) {
    has_coords <- tryCatch({
      graph$nodes[[0]]$get("coords", NULL) != reticulate::py_none()
    }, error = function(e) FALSE)
    if (isTRUE(has_coords)) {
      n <- as.integer(graph$number_of_nodes())
      coords_list <- lapply(0:(n - 1L), function(i) graph$nodes[[i]]$get("coords"))
      coords <- do.call(rbind, lapply(coords_list, function(v) as.numeric(v)))
    }
  }

  new_graph(A, coords = coords, normalized = normalized, directed = directed)
}

#' Export a graph to graph-tool
#'
#' Converts an `rgsp` graph to a Python graph-tool graph via `reticulate`.
#' Edge weights are stored under edge property `"weight"`.
#'
#' @param g A `gsp_graph`.
#' @return A graph-tool graph object (Python).
#' @export
graph_to_graphtool <- function(g) {
  if (!inherits(g, "gsp_graph")) stop("'g' must be a gsp_graph", call. = FALSE)
  .require_py_module("graph_tool")

  gt <- reticulate::import("graph_tool", delay_load = FALSE)
  graph <- gt$Graph(directed = isTRUE(g$directed))
  graph$add_vertex(as.integer(g$n))

  edges <- if (isTRUE(g$directed)) {
    Matrix::summary(g$adjacency)
  } else {
    Matrix::summary(Matrix::triu(g$adjacency, 1))
  }

  if (nrow(edges) > 0) {
    graph$add_edge_list(reticulate::r_to_py(cbind(edges$i - 1L, edges$j - 1L), convert = TRUE))
    wprop <- graph$new_edge_property("double")
    wprop$get_array()[seq_len(nrow(edges))] <- as.numeric(edges$x)
    graph$edge_properties[["weight"]] <- wprop
  }

  graph
}

#' Import a graph from graph-tool
#'
#' Converts a Python graph-tool graph to an `rgsp` `gsp_graph` via `reticulate`.
#'
#' @param graph A graph-tool graph object (Python).
#' @param weight Edge property name to use as weights (default `"weight"`).
#' @param normalized Logical; whether the returned graph uses a normalized Laplacian.
#' @return A `gsp_graph`.
#' @export
graph_from_graphtool <- function(graph, weight = "weight", normalized = FALSE) {
  .require_py_module("graph_tool")

  n <- as.integer(graph$num_vertices())
  directed <- isTRUE(graph$is_directed())

  # edge list
  edges <- reticulate::py_to_r(graph$get_edges())
  if (is.null(edges) || nrow(edges) == 0) {
    A <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    return(new_graph(A, normalized = normalized, directed = directed))
  }

  i <- as.integer(edges[, 1]) + 1L
  j <- as.integer(edges[, 2]) + 1L

  w <- NULL
  if (!is.null(weight) && !is.null(graph$edge_properties[[weight]])) {
    w <- as.numeric(reticulate::py_to_r(graph$edge_properties[[weight]]$get_array()))
  } else {
    w <- rep(1, length(i))
  }

  A <- Matrix::sparseMatrix(i = i, j = j, x = w, dims = c(n, n))
  if (!directed) {
    A <- Matrix::forceSymmetric(A, uplo = "U")
  }
  new_graph(Matrix::drop0(A), normalized = normalized, directed = directed)
}

#' Save a graph using NetworkX / graph-tool
#'
#' Convenience wrapper around Python NetworkX / graph-tool graph writers.
#'
#' @param g A `gsp_graph`.
#' @param path Output path.
#' @param fmt Format: `"graphml"`, `"gexf"`, or `"gml"`.
#' @param backend Backend: `"networkx"` or `"graph-tool"` (default `"networkx"`).
#' @param include_coords If `TRUE`, include `coords` as node attributes (may not
#'   be supported by all formats/writers).
#' @export
graph_save <- function(g,
                       path,
                       fmt = c("graphml", "gexf", "gml"),
                       backend = c("networkx", "graph-tool"),
                       include_coords = FALSE) {
  fmt <- match.arg(fmt)
  backend <- match.arg(backend)

  if (backend == "networkx") {
    .require_py_module("networkx")
    nx <- reticulate::import("networkx", delay_load = FALSE)
    g_nx <- graph_to_networkx(g, include_coords = include_coords)
    writer <- switch(
      fmt,
      graphml = nx$write_graphml,
      gexf = nx$write_gexf,
      gml = nx$write_gml
    )
    writer(g_nx, path)
    return(invisible(path))
  }

  .require_py_module("graph_tool")
  gt <- reticulate::import("graph_tool", delay_load = FALSE)
  g_gt <- graph_to_graphtool(g)
  g_gt$save(path, fmt = fmt)
  invisible(path)
}

#' Load a graph using NetworkX / graph-tool
#'
#' Convenience wrapper around Python NetworkX / graph-tool readers.
#'
#' @param path Path to input graph file.
#' @param fmt Format: `"graphml"`, `"gexf"`, or `"gml"`.
#' @param backend Backend: `"networkx"` or `"graph-tool"` (default `"networkx"`).
#' @param weight Weight attribute/property name (default `"weight"`).
#' @param normalized Logical; whether the returned graph uses a normalized Laplacian.
#' @export
graph_load <- function(path,
                       fmt = c("graphml", "gexf", "gml"),
                       backend = c("networkx", "graph-tool"),
                       weight = "weight",
                       normalized = FALSE) {
  fmt <- match.arg(fmt)
  backend <- match.arg(backend)

  if (backend == "networkx") {
    .require_py_module("networkx")
    nx <- reticulate::import("networkx", delay_load = FALSE)
    reader <- switch(
      fmt,
      graphml = nx$read_graphml,
      gexf = nx$read_gexf,
      gml = nx$read_gml
    )
    g_nx <- reader(path)
    return(graph_from_networkx(g_nx, weight = weight, normalized = normalized))
  }

  .require_py_module("graph_tool")
  gt <- reticulate::import("graph_tool", delay_load = FALSE)
  g_gt <- gt$load_graph(path, fmt = fmt)
  graph_from_graphtool(g_gt, weight = weight, normalized = normalized)
}

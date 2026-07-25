.cell_counts <- function(counts) {
  if (inherits(counts, "cudasparse")) {
    counts <- cudasparsr::to_dgCMatrix(counts)
  } else if (is.matrix(counts)) {
    if (!is.numeric(counts)) {
      stop("`counts` must contain numeric values.", call. = FALSE)
    }
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  } else if (!methods::is(counts, "Matrix")) {
    stop("`counts` must be a matrix, Matrix object, or cudasparse matrix.",
         call. = FALSE)
  }
  counts <- methods::as(
    methods::as(methods::as(counts, "dMatrix"), "generalMatrix"),
    "CsparseMatrix"
  )
  if (nrow(counts) < 2L || ncol(counts) < 2L) {
    stop("`counts` must contain at least two features and two cells.",
         call. = FALSE)
  }
  if (anyNA(counts@x) || any(!is.finite(counts@x)) || any(counts@x < 0)) {
    stop("`counts` must contain finite, non-negative values.",
         call. = FALSE)
  }
  counts
}

#' Normalize a single-cell count matrix
#'
#' Counts are expected with features in rows and cells in columns.
#'
#' @param counts A base matrix, `Matrix` sparse matrix, or `cudasparse`.
#' @param scale_factor Target library size.
#' @param log1p Whether to apply `log1p()` to non-zero normalized values.
#' @return A sparse `Matrix::dgCMatrix`.
#' @export
#' @examples
#' counts <- Matrix::Matrix(matrix(c(1, 0, 3, 2, 4, 1), 3), sparse = TRUE)
#' cuda_normalize_counts(counts)
cuda_normalize_counts <- function(counts, scale_factor = 10000,
                                  log1p = TRUE) {
  counts <- .cell_counts(counts)
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      is.na(scale_factor) || !is.finite(scale_factor) ||
      scale_factor <= 0) {
    stop("`scale_factor` must be a positive finite number.", call. = FALSE)
  }
  library_size <- as.numeric(Matrix::colSums(counts))
  if (any(library_size == 0)) {
    stop("Every cell must have a positive library size.", call. = FALSE)
  }
  normalized <- counts %*% Matrix::Diagonal(
    x = scale_factor / library_size
  )
  normalized <- methods::as(normalized, "dgCMatrix")
  if (isTRUE(log1p)) {
    normalized@x <- base::log1p(normalized@x)
  }
  normalized
}

#' Select highly variable features
#'
#' @param counts Feature-by-cell count or normalized expression matrix.
#' @param n_top Number of features to select.
#' @param min_mean Minimum feature mean.
#' @return A data frame ordered by dispersion with feature statistics and a
#'   `selected` indicator.
#' @export
#' @examples
#' set.seed(1)
#' counts <- Matrix::rsparsematrix(20, 10, density = 0.3)
#' counts@x <- abs(round(counts@x * 5))
#' cuda_hvg(counts, n_top = 5)
cuda_hvg <- function(counts, n_top = 2000L, min_mean = 0) {
  counts <- .cell_counts(counts)
  if (!is.numeric(n_top) || length(n_top) != 1L || is.na(n_top) ||
      n_top < 1 || n_top > nrow(counts) ||
      n_top != as.integer(n_top)) {
    stop("`n_top` must be a whole number between 1 and nrow(counts).",
         call. = FALSE)
  }
  if (!is.numeric(min_mean) || length(min_mean) != 1L ||
      is.na(min_mean) || !is.finite(min_mean) || min_mean < 0) {
    stop("`min_mean` must be a finite non-negative number.",
         call. = FALSE)
  }

  n_cells <- ncol(counts)
  mean_expression <- as.numeric(Matrix::rowSums(counts)) / n_cells
  squared <- counts
  squared@x <- squared@x^2
  mean_square <- as.numeric(Matrix::rowSums(squared)) / n_cells
  variance <- pmax(
    (mean_square - mean_expression^2) * n_cells / (n_cells - 1),
    0
  )
  dispersion <- variance / pmax(mean_expression, .Machine$double.eps)
  feature <- rownames(counts)
  if (is.null(feature)) {
    feature <- paste0("feature_", seq_len(nrow(counts)))
  }
  result <- data.frame(
    feature = feature,
    mean = mean_expression,
    variance = variance,
    dispersion = dispersion,
    selected = FALSE,
    stringsAsFactors = FALSE
  )
  eligible <- which(result$mean >= min_mean)
  ranked <- eligible[order(result$dispersion[eligible], decreasing = TRUE)]
  selected <- utils::head(ranked, as.integer(n_top))
  result$selected[selected] <- TRUE
  result[order(result$dispersion, decreasing = TRUE), , drop = FALSE]
}

.cell_n_components <- function(n_components) {
  if (!is.numeric(n_components) || length(n_components) != 1L ||
      is.na(n_components) || !is.finite(n_components) ||
      n_components < 1 || n_components != as.integer(n_components)) {
    stop("`n_components` must be a positive whole number.", call. = FALSE)
  }
  as.integer(n_components)
}

.cell_pca_from_normalized <- function(normalized, variable, n_components,
                                      scale., device) {
  n_components <- .cell_n_components(n_components)
  features <- variable$feature[variable$selected]
  feature_index <- match(features, rownames(normalized))
  if (anyNA(feature_index)) {
    feature_index <- match(features, paste0("feature_", seq_len(nrow(normalized))))
  }
  dense <- t(as.matrix(normalized[feature_index, , drop = FALSE]))
  max_components <- min(nrow(dense) - 1L, ncol(dense))
  if (n_components > max_components) {
    stop(
      sprintf("`n_components` cannot exceed %s for these data.",
              max_components),
      call. = FALSE
    )
  }
  fit <- cudalearnr::cuda_pca(
    dense,
    n_components = n_components,
    scale. = scale.,
    device = device
  )
  fit$features <- features
  fit
}

#' Run PCA on highly variable single-cell features
#'
#' @param counts Feature-by-cell matrix.
#' @param n_components Number of principal components.
#' @param n_hvg Number of highly variable features.
#' @param scale. Whether to scale selected features before PCA.
#' @param device Device passed to [cudalearnr::cuda_pca()].
#' @return A `cuda_pca` object with an additional `features` element.
#' @export
#' @examples
#' set.seed(2)
#' counts <- matrix(rpois(30 * 20, lambda = 2), 30, 20)
#' cuda_cell_pca(counts, n_components = 3, n_hvg = 10, device = "cpu")
cuda_cell_pca <- function(counts, n_components = 30L, n_hvg = 2000L,
                          scale. = TRUE,
                          device = c("auto", "cuda", "cpu")) {
  counts <- .cell_counts(counts)
  if (!is.numeric(n_hvg) || length(n_hvg) != 1L || is.na(n_hvg) ||
      !is.finite(n_hvg) || n_hvg < 1 || n_hvg != as.integer(n_hvg)) {
    stop("`n_hvg` must be a positive whole number.", call. = FALSE)
  }
  n_hvg <- min(as.integer(n_hvg), nrow(counts))
  n_components <- .cell_n_components(n_components)
  normalized <- cuda_normalize_counts(counts)
  variable <- cuda_hvg(normalized, n_top = n_hvg)
  .cell_pca_from_normalized(
    normalized = normalized,
    variable = variable,
    n_components = n_components,
    scale. = scale.,
    device = device
  )
}

#' Construct a nearest-neighbour representation
#'
#' @param embedding Cell-by-component matrix or `cuda_pca` result.
#' @param k Number of neighbours.
#' @param metric Distance metric.
#' @param device Computation device.
#' @param batch_size Maximum query rows per exact kNN distance block. Smaller
#'   values reduce peak distance memory.
#' @return A `cuda_knn` object.
#' @export
#' @examples
#' cuda_cell_neighbors(
#'   matrix(rnorm(60), 20, 3),
#'   k = 5,
#'   batch_size = 4,
#'   device = "cpu"
#' )
cuda_cell_neighbors <- function(embedding, k = 15L,
                                metric = c("euclidean", "cosine"),
                                device = c("auto", "cuda", "cpu"),
                                batch_size = 256L) {
  if (inherits(embedding, "cuda_pca")) {
    embedding <- embedding$x
  }
  cudalearnr::cuda_knn(
    embedding,
    k = k,
    metric = match.arg(metric),
    device = device,
    batch_size = batch_size
  )
}

#' Run the initial cudacellr workflow
#'
#' @param counts Feature-by-cell count matrix.
#' @param n_hvg Number of variable features.
#' @param n_components PCA components.
#' @param k Neighbours per cell.
#' @param device Computation device.
#' @param batch_size Maximum query rows per exact kNN distance block.
#' @return A list containing normalized counts, variable features, PCA, and kNN.
#' @export
#' @examples
#' set.seed(3)
#' counts <- matrix(rpois(40 * 25, 2), 40, 25)
#' cudacell_workflow(
#'   counts, n_hvg = 15, n_components = 5, k = 5, device = "cpu"
#' )
cudacell_workflow <- function(counts, n_hvg = 2000L,
                              n_components = 30L, k = 15L,
                              device = c("auto", "cuda", "cpu"),
                              batch_size = 256L) {
  counts <- .cell_counts(counts)
  if (!is.numeric(n_hvg) || length(n_hvg) != 1L || is.na(n_hvg) ||
      !is.finite(n_hvg) || n_hvg < 1 || n_hvg != as.integer(n_hvg)) {
    stop("`n_hvg` must be a positive whole number.", call. = FALSE)
  }
  n_components <- .cell_n_components(n_components)
  normalized <- cuda_normalize_counts(counts)
  variable <- cuda_hvg(
    normalized,
    n_top = min(as.integer(n_hvg), nrow(counts))
  )
  pca <- .cell_pca_from_normalized(
    normalized = normalized,
    variable = variable,
    n_components = n_components,
    scale. = TRUE,
    device = device
  )
  neighbours <- cuda_cell_neighbors(
    pca,
    k = k,
    device = device,
    batch_size = batch_size
  )
  structure(
    list(
      normalized = normalized,
      variable_features = variable,
      pca = pca,
      neighbors = neighbours
    ),
    class = "cudacell_workflow"
  )
}

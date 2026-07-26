.cell_provenance_metadata <- function(stages) {
  provenance <- cudatensr::cuda_provenance(stages)
  list(
    provenance_schema = attr(provenance, "schema", exact = TRUE),
    compute_device = attr(provenance, "compute_device", exact = TRUE),
    compute_stages = attr(provenance, "compute_stages", exact = TRUE)
  )
}

.with_cell_provenance <- function(x, stages) {
  metadata <- .cell_provenance_metadata(stages)
  if (is.list(x) && is.null(dim(x)) && !methods::is(x, "Matrix")) {
    x$provenance_schema <- metadata$provenance_schema
    x$compute_device <- metadata$compute_device
    x$compute_stages <- metadata$compute_stages
    return(x)
  }
  attr(x, "provenance_schema") <- metadata$provenance_schema
  attr(x, "compute_device") <- metadata$compute_device
  attr(x, "compute_stages") <- metadata$compute_stages
  x
}

.cell_has_compute_stages <- function(x) {
  if (is.list(x) && "compute_stages" %in% names(x)) {
    return(TRUE)
  }
  !is.null(attr(x, "compute_stages", exact = TRUE))
}

.cell_source_stages <- function(x) {
  if (!.cell_has_compute_stages(x)) {
    return(list())
  }
  attr(
    cudatensr::cuda_provenance(x),
    "compute_stages",
    exact = TRUE
  )
}

.cell_cpu_stage <- function(backend = "Matrix",
                            reason = "algorithm_cpu_only") {
  cudatensr::cuda_stage(
    requested_device = "fixed-cpu",
    device = "cpu",
    backend = backend,
    selection_reason = reason,
    fallback = FALSE,
    output_device = "cpu"
  )
}

.cell_prefix_stages <- function(stages, prefix) {
  if (!length(stages)) {
    return(stages)
  }
  names(stages) <- paste0(prefix, names(stages))
  stages
}

.cell_append_stages <- function(stages, additions) {
  duplicate <- intersect(names(stages), names(additions))
  if (length(duplicate)) {
    stop(
      sprintf(
        "Duplicate compute stage name: %s.",
        paste(duplicate, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  c(stages, additions)
}

.cell_input_stages <- function(x) {
  stages <- .cell_source_stages(x)
  if (inherits(x, "cudasparse")) {
    stages <- .cell_prefix_stages(stages, "source_")
  }
  if (inherits(x, "cudasparse") && identical(x$device, "cuda")) {
    stages <- .cell_append_stages(
      stages,
      list(
        input_materialization = .cell_cpu_stage(
          reason = "input_transfer"
        )
      )
    )
  }
  stages
}

#' @importFrom cudatensr cuda_provenance
#' @export
cudatensr::cuda_provenance

#' Inspect provenance stored on native single-cell containers
#'
#' These methods expose the shared [cudatensr::cuda_provenance()] contract for
#' `SingleCellExperiment` and `Seurat` objects produced by `cudacellr`.
#'
#' @param x A supported native single-cell container.
#' @return A `cuda_provenance` data frame.
#' @name cuda_provenance-methods
NULL

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
  count_dimnames <- dimnames(counts)
  counts <- methods::as(
    methods::as(methods::as(counts, "dMatrix"), "generalMatrix"),
    "CsparseMatrix"
  )
  dimnames(counts) <- count_dimnames
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
#' @return A sparse `Matrix::dgCMatrix` with the input feature and cell names.
#' @export
#' @examples
#' counts <- Matrix::Matrix(matrix(c(1, 0, 3, 2, 4, 1), 3), sparse = TRUE)
#' cuda_normalize_counts(counts)
cuda_normalize_counts <- function(counts, scale_factor = 10000,
                                  log1p = TRUE) {
  stages <- .cell_input_stages(counts)
  counts <- .cell_counts(counts)
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      is.na(scale_factor) || !is.finite(scale_factor) ||
      scale_factor <= 0) {
    stop("`scale_factor` must be a positive finite number.", call. = FALSE)
  }
  if (!is.logical(log1p) || length(log1p) != 1L || is.na(log1p)) {
    stop("`log1p` must be TRUE or FALSE.", call. = FALSE)
  }
  library_size <- as.numeric(Matrix::colSums(counts))
  if (any(library_size == 0)) {
    stop("Every cell must have a positive library size.", call. = FALSE)
  }
  normalized <- counts %*% Matrix::Diagonal(
    x = scale_factor / library_size
  )
  normalized <- methods::as(normalized, "dgCMatrix")
  dimnames(normalized) <- dimnames(counts)
  if (log1p) {
    normalized@x <- base::log1p(normalized@x)
  }
  stages <- .cell_append_stages(
    stages,
    list(normalization = .cell_cpu_stage())
  )
  .with_cell_provenance(normalized, stages)
}

#' Select highly variable features
#'
#' @param counts Feature-by-cell count or normalized expression matrix.
#' @param n_top Number of features to select.
#' @param min_mean Minimum feature mean.
#' @return A data frame ordered by dispersion with feature statistics and a
#'   `selected` indicator. The `index` column retains the original feature
#'   position even when feature names are duplicated.
#' @export
#' @examples
#' set.seed(1)
#' counts <- Matrix::rsparsematrix(20, 10, density = 0.3)
#' counts@x <- abs(round(counts@x * 5))
#' cuda_hvg(counts, n_top = 5)
cuda_hvg <- function(counts, n_top = 2000L, min_mean = 0) {
  stages <- .cell_input_stages(counts)
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
    index = seq_len(nrow(counts)),
    stringsAsFactors = FALSE
  )
  eligible <- which(result$mean >= min_mean)
  ranked <- eligible[order(result$dispersion[eligible], decreasing = TRUE)]
  selected <- utils::head(ranked, as.integer(n_top))
  result$selected[selected] <- TRUE
  result <- result[order(result$dispersion, decreasing = TRUE), , drop = FALSE]
  stages <- .cell_append_stages(
    stages,
    list(hvg = .cell_cpu_stage())
  )
  .with_cell_provenance(result, stages)
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
  selected <- variable[variable$selected, , drop = FALSE]
  features <- selected$feature
  feature_index <- selected$index
  dense <- t(as.matrix(normalized[feature_index, , drop = FALSE]))
  rownames(dense) <- colnames(normalized)
  colnames(dense) <- features
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
  preprocessing <- .cell_source_stages(variable)
  pca_stages <- .cell_prefix_stages(
    .cell_source_stages(fit),
    "pca_"
  )
  .with_cell_provenance(
    fit,
    .cell_append_stages(preprocessing, pca_stages)
  )
}

#' Run PCA on highly variable single-cell features
#'
#' @param counts Feature-by-cell matrix.
#' @param n_components Number of principal components.
#' @param n_hvg Number of highly variable features.
#' @param scale. Whether to scale selected features before PCA.
#' @param device Device passed to [cudalearnr::cuda_pca()].
#' @param scale_factor Target library size used during normalization.
#' @param log1p Whether to apply `log1p()` after library-size normalization.
#' @param min_mean Minimum feature mean used during highly variable feature
#'   selection.
#' @return A `cuda_pca` object with an additional `features` element. Cell
#'   names label PCA score rows and selected feature names label loadings.
#' @export
#' @examples
#' set.seed(2)
#' counts <- matrix(rpois(30 * 20, lambda = 2), 30, 20)
#' cuda_cell_pca(counts, n_components = 3, n_hvg = 10, device = "cpu")
cuda_cell_pca <- function(counts, n_components = 30L, n_hvg = 2000L,
                          scale. = TRUE,
                          device = c("auto", "cuda", "cpu"),
                          scale_factor = 10000, log1p = TRUE,
                          min_mean = 0) {
  device <- match.arg(device)
  cudatensr::cuda_select_device(device)
  input_stages <- .cell_input_stages(counts)
  counts <- .cell_counts(counts)
  if (!is.numeric(n_hvg) || length(n_hvg) != 1L || is.na(n_hvg) ||
      !is.finite(n_hvg) || n_hvg < 1 || n_hvg != as.integer(n_hvg)) {
    stop("`n_hvg` must be a positive whole number.", call. = FALSE)
  }
  n_hvg <- min(as.integer(n_hvg), nrow(counts))
  n_components <- .cell_n_components(n_components)
  normalized <- cuda_normalize_counts(
    counts,
    scale_factor = scale_factor,
    log1p = log1p
  )
  if (length(input_stages)) {
    normalized <- .with_cell_provenance(
      normalized,
      .cell_append_stages(
        input_stages,
        .cell_source_stages(normalized)
      )
    )
  }
  variable <- cuda_hvg(
    normalized,
    n_top = n_hvg,
    min_mean = min_mean
  )
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
#' @return A `cuda_knn` object retaining cell names from the embedding.
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
  device <- match.arg(device)
  cudatensr::cuda_select_device(device)
  source_stages <- .cell_source_stages(embedding)
  if (inherits(embedding, "cuda_pca")) {
    embedding <- embedding$x
  }
  result <- cudalearnr::cuda_knn(
    embedding,
    k = k,
    metric = match.arg(metric),
    device = device,
    batch_size = batch_size
  )
  neighbor_stages <- .cell_prefix_stages(
    .cell_source_stages(result),
    "knn_"
  )
  .with_cell_provenance(
    result,
    .cell_append_stages(source_stages, neighbor_stages)
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
#' @param scale_factor Target library size used during normalization.
#' @param log1p Whether to apply `log1p()` after library-size normalization.
#' @param min_mean Minimum feature mean used during highly variable feature
#'   selection.
#' @param scale. Whether to scale selected features before PCA.
#' @return A list containing normalized counts, variable features, PCA, and
#'   kNN. Feature and cell identifiers are retained throughout all stages.
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
                              batch_size = 256L,
                              scale_factor = 10000, log1p = TRUE,
                              min_mean = 0, scale. = TRUE) {
  device <- match.arg(device)
  cudatensr::cuda_select_device(device)
  input_stages <- .cell_input_stages(counts)
  counts <- .cell_counts(counts)
  if (!is.numeric(n_hvg) || length(n_hvg) != 1L || is.na(n_hvg) ||
      !is.finite(n_hvg) || n_hvg < 1 || n_hvg != as.integer(n_hvg)) {
    stop("`n_hvg` must be a positive whole number.", call. = FALSE)
  }
  n_components <- .cell_n_components(n_components)
  normalized <- cuda_normalize_counts(
    counts,
    scale_factor = scale_factor,
    log1p = log1p
  )
  if (length(input_stages)) {
    normalized <- .with_cell_provenance(
      normalized,
      .cell_append_stages(
        input_stages,
        .cell_source_stages(normalized)
      )
    )
  }
  variable <- cuda_hvg(
    normalized,
    n_top = min(as.integer(n_hvg), nrow(counts)),
    min_mean = min_mean
  )
  pca <- .cell_pca_from_normalized(
    normalized = normalized,
    variable = variable,
    n_components = n_components,
    scale. = scale.,
    device = device
  )
  neighbours <- cuda_cell_neighbors(
    pca,
    k = k,
    device = device,
    batch_size = batch_size
  )
  output <- structure(
    list(
      normalized = normalized,
      variable_features = variable,
      pca = pca,
      neighbors = neighbours
    ),
    class = "cudacell_workflow"
  )
  .with_cell_provenance(output, .cell_source_stages(neighbours))
}

#' @export
print.cudacell_workflow <- function(x, ...) {
  cat(sprintf(
    paste0(
      "<cudacell_workflow features=%s cells=%s hvg=%s ",
      "components=%s k=%s pca_device=%s compute=%s>\n"
    ),
    nrow(x$normalized),
    ncol(x$normalized),
    sum(x$variable_features$selected),
    ncol(x$pca$x),
    ncol(x$neighbors$index),
    x$pca$device,
    x$compute_device
  ))
  invisible(x)
}

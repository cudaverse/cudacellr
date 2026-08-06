.sce_schema <- "cudacellr-sce/1"

.sce_stop <- function(message, class, ...) {
  fields <- list(...)
  condition <- structure(
    c(
      list(message = message, call = NULL),
      fields
    ),
    class = c(class, "cudacell_sce_error", "error", "condition")
  )
  stop(condition)
}

.sce_require <- function() {
  packages <- c(
    "SingleCellExperiment",
    "SummarizedExperiment",
    "S4Vectors"
  )
  available <- vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  if (!all(available)) {
    .sce_stop(
      paste0(
        "Install the Bioconductor package 'SingleCellExperiment' to use ",
        "SingleCellExperiment adapters. Run ",
        "`BiocManager::install(\"SingleCellExperiment\")`."
      ),
      class = "cudacell_sce_dependency_error",
      missing_packages = packages[!available]
    )
  }
  invisible(TRUE)
}

.sce_scalar_name <- function(x, argument, allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !nzchar(x) || !identical(trimws(x), x)) {
    .sce_stop(
      sprintf("`%s` must be one non-empty character string.", argument),
      class = "cudacell_sce_name_error",
      argument = argument
    )
  }
  x
}

.sce_flag <- function(x, argument) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .sce_stop(
      sprintf("`%s` must be TRUE or FALSE.", argument),
      class = "cudacell_sce_name_error",
      argument = argument
    )
  }
  x
}

.sce_row_fields <- function(prefix) {
  paste0(
    prefix,
    "_",
    c("mean", "variance", "dispersion", "hvg", "hvg_rank")
  )
}

.sce_size_factor_field <- function(prefix) {
  paste0(prefix, "_size_factor")
}

.sce_conflicts <- function(x, normalized_assay, reduced_dim, col_pair,
                           row_data_prefix, set_size_factors) {
  conflicts <- character()
  if (normalized_assay %in% SummarizedExperiment::assayNames(x)) {
    conflicts <- c(
      conflicts,
      sprintf("assay '%s'", normalized_assay)
    )
  }
  if (reduced_dim %in% SingleCellExperiment::reducedDimNames(x)) {
    conflicts <- c(
      conflicts,
      sprintf("reducedDim '%s'", reduced_dim)
    )
  }
  if (col_pair %in% SingleCellExperiment::colPairNames(x)) {
    conflicts <- c(
      conflicts,
      sprintf("colPair '%s'", col_pair)
    )
  }
  row_conflicts <- intersect(
    .sce_row_fields(row_data_prefix),
    names(SummarizedExperiment::rowData(x))
  )
  if (length(row_conflicts)) {
    conflicts <- c(
      conflicts,
      sprintf("rowData '%s'", row_conflicts)
    )
  }
  size_factor_field <- .sce_size_factor_field(row_data_prefix)
  if (size_factor_field %in% names(SummarizedExperiment::colData(x))) {
    conflicts <- c(
      conflicts,
      sprintf("colData '%s'", size_factor_field)
    )
  }
  if ("cudacellr" %in% names(S4Vectors::metadata(x))) {
    conflicts <- c(conflicts, "metadata 'cudacellr'")
  }
  if (set_size_factors &&
      !is.null(SingleCellExperiment::sizeFactors(x))) {
    conflicts <- c(conflicts, "canonical sizeFactors")
  }
  conflicts
}

.sce_assay <- function(x, assay, realize) {
  assay_names <- SummarizedExperiment::assayNames(x)
  if (!assay %in% assay_names) {
    available <- if (length(assay_names)) {
      paste(sprintf("'%s'", assay_names), collapse = ", ")
    } else {
      "<none>"
    }
    .sce_stop(
      sprintf(
        "`assay` '%s' was not found. Available assays: %s.",
        assay,
        available
      ),
      class = "cudacell_sce_assay_error",
      assay = assay,
      available_assays = assay_names
    )
  }
  counts <- SummarizedExperiment::assay(
    x,
    assay,
    withDimnames = TRUE
  )
  assay_class <- class(counts)
  materialized <- FALSE
  if (inherits(counts, "DelayedMatrix")) {
    if (!realize) {
      .sce_stop(
        paste0(
          "The selected assay is a DelayedMatrix. cudacellr will not ",
          "silently realize it in memory; set `realize = TRUE` after ",
          "confirming that the assay fits in memory."
        ),
        class = "cudacell_sce_assay_error",
        assay = assay,
        assay_class = class(counts)
      )
    }
    original_dimnames <- dimnames(counts)
    counts <- tryCatch(
      methods::as(counts, "CsparseMatrix"),
      error = function(error) {
        .sce_stop(
          sprintf(
            "Could not realize assay '%s' as a sparse Matrix: %s",
            assay,
            conditionMessage(error)
          ),
          class = "cudacell_sce_assay_error",
          assay = assay,
          assay_class = class(counts)
        )
      }
    )
    dimnames(counts) <- original_dimnames
    materialized <- TRUE
  }
  list(
    counts = counts,
    class = assay_class,
    materialized = materialized
  )
}

.sce_hvg_columns <- function(variable_features, n_features) {
  required <- c(
    "index",
    "mean",
    "variance",
    "dispersion",
    "selected"
  )
  if (!is.data.frame(variable_features) ||
      !all(required %in% names(variable_features))) {
    .sce_stop(
      "The workflow returned incomplete feature statistics.",
      class = "cudacell_sce_mapping_error"
    )
  }
  positions <- variable_features$index
  if (length(positions) != n_features ||
      !is.numeric(positions) ||
      anyNA(positions) ||
      any(!is.finite(positions)) ||
      any(positions != as.integer(positions)) ||
      !identical(sort(as.integer(positions)), seq_len(n_features))) {
    .sce_stop(
      "The workflow returned an invalid feature-position mapping.",
      class = "cudacell_sce_mapping_error"
    )
  }
  output <- list(
    mean = numeric(n_features),
    variance = numeric(n_features),
    dispersion = numeric(n_features),
    hvg = logical(n_features),
    hvg_rank = integer(n_features)
  )
  output$mean[positions] <- variable_features$mean
  output$variance[positions] <- variable_features$variance
  output$dispersion[positions] <- variable_features$dispersion
  output$hvg[positions] <- variable_features$selected
  output$hvg_rank[positions] <- seq_len(n_features)
  output
}

.sce_col_pair <- function(neighbors, n_cells) {
  index <- neighbors$index
  distance <- neighbors$distance
  if (!is.matrix(index) || !is.numeric(index) ||
      !is.matrix(distance) || !is.numeric(distance) ||
      !identical(dim(index), dim(distance)) ||
      nrow(index) != n_cells) {
    .sce_stop(
      "The workflow returned an invalid nearest-neighbour mapping.",
      class = "cudacell_sce_mapping_error"
    )
  }
  k <- ncol(index)
  from <- rep(seq_len(n_cells), each = k)
  raw_to <- as.vector(t(index))
  distances <- as.numeric(as.vector(t(distance)))
  ranks <- rep.int(seq_len(k), times = n_cells)
  if (anyNA(raw_to) || any(!is.finite(raw_to)) ||
      any(raw_to != as.integer(raw_to))) {
    .sce_stop(
      "The workflow returned non-integer neighbour indices.",
      class = "cudacell_sce_mapping_error"
    )
  }
  to <- as.integer(raw_to)
  if (any(to < 1L | to > n_cells) ||
      anyNA(distances) || any(!is.finite(distances)) ||
      any(distances < 0) ||
      any(from == to) ||
      anyDuplicated(paste(from, to, sep = ":"))) {
    .sce_stop(
      "The workflow returned invalid neighbour indices or distances.",
      class = "cudacell_sce_mapping_error"
    )
  }
  S4Vectors::SelfHits(
    from = from,
    to = to,
    nnode = n_cells,
    distance = distances,
    rank = ranks
  )
}

.sce_provenance <- function(x) {
  .sce_require()
  record <- S4Vectors::metadata(x)[["cudacellr"]]
  if (is.null(record)) {
    .sce_stop(
      paste0(
        "This SingleCellExperiment has no cudacellr metadata. ",
        "Run `cudacell_sce()` first."
      ),
      class = "cudacell_sce_metadata_error"
    )
  }
  if (!is.list(record) ||
      !identical(record$schema, .sce_schema) ||
      is.null(record$provenance_schema) ||
      is.null(record$compute_device) ||
      is.null(record$compute_stages)) {
    .sce_stop(
      "The stored cudacellr SingleCellExperiment metadata is invalid.",
      class = "cudacell_sce_metadata_error"
    )
  }
  holder <- list(
    provenance_schema = record$provenance_schema,
    compute_device = record$compute_device,
    compute_stages = record$compute_stages
  )
  cudaverse::cuda_provenance(holder)
}

#' @rdname cuda_provenance-methods
#' @method cuda_provenance SingleCellExperiment
#' @export
cuda_provenance.SingleCellExperiment <- function(x) {
  .sce_provenance(x)
}

#' Run a cudacell workflow on a SingleCellExperiment
#'
#' `cudacell_sce()` reads one assay, runs [cudacell_workflow()], and returns a
#' modified copy of the input object. Existing assays, row and column
#' metadata, reduced dimensions, alternative experiments, pairings, and
#' top-level metadata are retained. New results are stored in native
#' `SingleCellExperiment` locations:
#'
#' - normalized expression in an assay;
#' - PCA scores in a reduced dimension;
#' - highly variable feature statistics in `rowData`;
#' - directed nearest-neighbour relationships in a `colPair`;
#' - parameters and compute provenance in `metadata(x)$cudacellr`.
#'
#' The selected assay must currently be a base matrix or `Matrix` object.
#' Delayed assays are never realized silently. Set `realize = TRUE` only after
#' confirming that the selected assay fits in memory.
#'
#' @param x A `SingleCellExperiment`.
#' @param assay Exact name of the feature-by-cell count assay.
#' @param n_hvg Number of highly variable features.
#' @param n_components Number of PCA components.
#' @param k Number of neighbours per cell.
#' @param device Computation device.
#' @param batch_size Maximum query rows per exact kNN distance block.
#' @param scale_factor Target library size used during normalization.
#' @param log1p Whether to apply `log1p()` after library-size normalization.
#' @param min_mean Minimum feature mean used for HVG selection.
#' @param scale. Whether to scale selected features before PCA.
#' @param normalized_assay Name for the added normalized assay. `NULL` chooses
#'   `"cudacell_logcounts"` when `log1p = TRUE` and
#'   `"cudacell_normalized"` otherwise.
#' @param reduced_dim Name for the added PCA reduced dimension.
#' @param col_pair Name for the added directed kNN `colPair`.
#' @param row_data_prefix Prefix for added feature statistics and the
#'   namespaced column size factor.
#' @param set_size_factors Whether to also replace the canonical
#'   [SingleCellExperiment::sizeFactors()] values. The namespaced
#'   `<row_data_prefix>_size_factor` column is always added to `colData`.
#' @param overwrite Whether explicitly named cudacellr output fields may be
#'   replaced. Other object content is never overwritten.
#' @param realize Whether a delayed assay may be explicitly materialized as a
#'   sparse `Matrix`.
#' @return A valid `SingleCellExperiment` of the same class as `x`.
#' @export
#' @examples
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   counts <- matrix(
#'     rpois(30 * 12, lambda = 2),
#'     nrow = 30,
#'     dimnames = list(paste0("gene_", 1:30), paste0("cell_", 1:12))
#'   )
#'   sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts)
#'   )
#'   sce <- cudacell_sce(
#'     sce,
#'     n_hvg = 10,
#'     n_components = 3,
#'     k = 3,
#'     device = "cpu"
#'   )
#'   SingleCellExperiment::reducedDim(sce, "CUDACELL_PCA")
#'   cuda_provenance(sce)
#' }
cudacell_sce <- function(x, assay = "counts", n_hvg = 2000L,
                         n_components = 30L, k = 15L,
                         device = c("auto", "cuda", "cpu"),
                         batch_size = 256L, scale_factor = 10000,
                         log1p = TRUE, min_mean = 0, scale. = TRUE,
                         normalized_assay = NULL,
                         reduced_dim = "CUDACELL_PCA",
                         col_pair = "CUDACELL_KNN",
                         row_data_prefix = "cudacell",
                         set_size_factors = FALSE,
                         overwrite = FALSE, realize = FALSE) {
  .sce_require()
  if (!methods::is(x, "SingleCellExperiment")) {
    .sce_stop(
      "`x` must be a SingleCellExperiment object.",
      class = "cudacell_sce_type_error"
    )
  }
  assay <- .sce_scalar_name(assay, "assay")
  log1p <- .sce_flag(log1p, "log1p")
  scale. <- .sce_flag(scale., "scale.")
  set_size_factors <- .sce_flag(
    set_size_factors,
    "set_size_factors"
  )
  overwrite <- .sce_flag(overwrite, "overwrite")
  realize <- .sce_flag(realize, "realize")
  if (is.null(normalized_assay)) {
    normalized_assay <- if (log1p) {
      "cudacell_logcounts"
    } else {
      "cudacell_normalized"
    }
  }
  normalized_assay <- .sce_scalar_name(
    normalized_assay,
    "normalized_assay"
  )
  if (identical(normalized_assay, assay)) {
    .sce_stop(
      "`normalized_assay` must differ from the input `assay`.",
      class = "cudacell_sce_name_error",
      argument = "normalized_assay"
    )
  }
  reduced_dim <- .sce_scalar_name(reduced_dim, "reduced_dim")
  col_pair <- .sce_scalar_name(col_pair, "col_pair")
  row_data_prefix <- .sce_scalar_name(
    row_data_prefix,
    "row_data_prefix"
  )

  assay_names <- SummarizedExperiment::assayNames(x)
  if (!assay %in% assay_names) {
    .sce_assay(x, assay, realize = FALSE)
  }
  conflicts <- .sce_conflicts(
    x,
    normalized_assay = normalized_assay,
    reduced_dim = reduced_dim,
    col_pair = col_pair,
    row_data_prefix = row_data_prefix,
    set_size_factors = set_size_factors
  )
  if (length(conflicts) && !overwrite) {
    .sce_stop(
      paste0(
        "cudacellr output fields already exist: ",
        paste(conflicts, collapse = ", "),
        ". Set `overwrite = TRUE` to replace only these fields."
      ),
      class = "cudacell_sce_collision_error",
      conflicts = conflicts
    )
  }

  device <- match.arg(device)
  cudaverse::cuda_select_device(device)
  input <- .sce_assay(x, assay, realize = realize)
  counts <- .cell_counts(input$counts)
  fit <- cudacell_workflow(
    counts,
    n_hvg = n_hvg,
    n_components = n_components,
    k = k,
    device = device,
    batch_size = batch_size,
    scale_factor = scale_factor,
    log1p = log1p,
    min_mean = min_mean,
    scale. = scale.
  )
  library_size <- as.numeric(Matrix::colSums(counts))
  size_factors <- library_size / scale_factor

  if (!identical(dim(fit$normalized), dim(x)) ||
      !identical(rownames(fit$normalized), rownames(x)) ||
      !identical(colnames(fit$normalized), colnames(x)) ||
      nrow(fit$pca$x) != ncol(x) ||
      !identical(rownames(fit$pca$x), colnames(x))) {
    .sce_stop(
      "The workflow results do not align with the SingleCellExperiment.",
      class = "cudacell_sce_mapping_error"
    )
  }
  hvg <- .sce_hvg_columns(fit$variable_features, nrow(x))
  neighbor_hits <- .sce_col_pair(fit$neighbors, ncol(x))

  out <- x
  SummarizedExperiment::assay(
    out,
    normalized_assay,
    withDimnames = TRUE
  ) <- fit$normalized

  row_data <- SummarizedExperiment::rowData(out)
  row_fields <- .sce_row_fields(row_data_prefix)
  for (index in seq_along(row_fields)) {
    row_data[[row_fields[[index]]]] <- hvg[[index]]
  }
  SummarizedExperiment::rowData(out) <- row_data

  col_data <- SummarizedExperiment::colData(out)
  size_factor_field <- .sce_size_factor_field(row_data_prefix)
  col_data[[size_factor_field]] <- size_factors
  SummarizedExperiment::colData(out) <- col_data

  SingleCellExperiment::reducedDim(
    out,
    reduced_dim,
    withDimnames = TRUE
  ) <- fit$pca$x
  SingleCellExperiment::colPair(out, col_pair) <- neighbor_hits
  if (set_size_factors) {
    SingleCellExperiment::sizeFactors(out) <- size_factors
  }

  object_metadata <- S4Vectors::metadata(out)
  object_metadata[["cudacellr"]] <- list(
    schema = .sce_schema,
    package_version = as.character(
      utils::packageVersion("cudacellr")
    ),
    input_assay = assay,
    input_assay_class = input$class,
    input_assay_materialized = input$materialized,
    normalized_assay = normalized_assay,
    reduced_dim = reduced_dim,
    col_pair = col_pair,
    row_data_fields = row_fields,
    size_factor_field = size_factor_field,
    canonical_size_factors = set_size_factors,
    parameters = list(
      n_hvg = as.integer(n_hvg),
      n_components = as.integer(n_components),
      k = as.integer(k),
      requested_device = device,
      batch_size = as.integer(batch_size),
      scale_factor = scale_factor,
      log1p = log1p,
      min_mean = min_mean,
      scale = scale.
    ),
    provenance_schema = fit$provenance_schema,
    compute_device = fit$compute_device,
    compute_stages = fit$compute_stages
  )
  S4Vectors::metadata(out) <- object_metadata

  methods::validObject(out)
  out
}

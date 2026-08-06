.skip_sce_dependencies <- function() {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
}

.sce_example <- function(with_names = TRUE) {
  .skip_sce_dependencies()
  counts <- matrix(
    c(
      4, 0, 1, 3, 2, 1, 0, 2, 1, 4, 2, 3,
      0, 2, 3, 1, 1, 4, 2, 0, 3, 1, 2, 4,
      5, 1, 0, 2, 3, 1, 4, 2, 1, 0, 3, 2,
      1, 4, 2, 1, 1, 3, 2, 4, 0, 2, 1, 3,
      2, 1, 4, 3, 1, 2, 0, 1, 4, 3, 2, 1,
      3, 2, 1, 4, 2, 1, 3, 0, 2, 4, 1, 2
    ),
    nrow = 6,
    byrow = TRUE
  )
  if (with_names) {
    rownames(counts) <- c(
      "duplicated_gene",
      "duplicated_gene",
      paste0("gene_", 3:6)
    )
    colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))
  }
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  second_assay <- counts
  second_assay@x <- second_assay@x + 10

  existing_embedding <- matrix(
    seq_len(ncol(counts) * 2L),
    nrow = ncol(counts),
    dimnames = list(colnames(counts), c("old_1", "old_2"))
  )
  attr(existing_embedding, "source") <- "preserve-me"
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      untouched = second_assay
    ),
    rowData = S4Vectors::DataFrame(
      symbol = paste0("symbol_", seq_len(nrow(counts))),
      category = factor(rep(c("a", "b"), length.out = nrow(counts)))
    ),
    colData = S4Vectors::DataFrame(
      batch = factor(rep(c("one", "two"), length.out = ncol(counts))),
      score = seq_len(ncol(counts))
    ),
    metadata = list(
      owner = "preserve-me",
      nested = list(value = 42L)
    ),
    reducedDims = list(existing = existing_embedding)
  )

  alternative <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      spike = matrix(
        seq_len(2L * ncol(counts)),
        nrow = 2L,
        dimnames = list(
          if (with_names) c("spike_1", "spike_2") else NULL,
          colnames(counts)
        )
      )
    ),
    rowData = S4Vectors::DataFrame(kind = c("a", "b")),
    metadata = list(source = "alternative")
  )
  SingleCellExperiment::altExp(sce, "spike") <- alternative
  SingleCellExperiment::mainExpName(sce) <- "endogenous"
  SingleCellExperiment::sizeFactors(sce) <- seq(
    0.75,
    1.25,
    length.out = ncol(sce)
  )
  SingleCellExperiment::colLabels(sce) <- factor(
    rep(c("A", "B"), length.out = ncol(sce))
  )
  SingleCellExperiment::rowPair(sce, "existing_rows") <-
    S4Vectors::SelfHits(1L, 2L, nnode = nrow(sce), value = 7)
  SingleCellExperiment::colPair(sce, "existing_cells") <-
    S4Vectors::SelfHits(1L, 2L, nnode = ncol(sce), value = 9)
  sce
}

.sce_workflow_args <- function() {
  list(
    n_hvg = 5,
    n_components = 3,
    k = 3,
    batch_size = 2,
    device = "cpu"
  )
}

test_that("cudacell_sce writes native results without losing SCE content", {
  sce <- .sce_example()
  original <- sce
  original_assays <- SummarizedExperiment::assays(
    sce,
    withDimnames = FALSE
  )
  original_row_data <- SummarizedExperiment::rowData(sce)
  original_col_data <- SummarizedExperiment::colData(sce)
  original_metadata <- S4Vectors::metadata(sce)
  original_reduced <- SingleCellExperiment::reducedDim(
    sce,
    "existing",
    withDimnames = FALSE
  )
  original_alt <- SingleCellExperiment::altExp(sce, "spike")
  original_row_pair <- SingleCellExperiment::rowPair(
    sce,
    "existing_rows"
  )
  original_col_pair <- SingleCellExperiment::colPair(
    sce,
    "existing_cells"
  )
  original_size_factors <- SingleCellExperiment::sizeFactors(sce)
  original_labels <- SingleCellExperiment::colLabels(sce)

  arguments <- .sce_workflow_args()
  out <- do.call(cudacell_sce, c(list(x = sce), arguments))
  fit <- do.call(
    cudacell_workflow,
    c(
      list(
        counts = SummarizedExperiment::assay(sce, "counts")
      ),
      arguments
    )
  )

  expect_identical(sce, original)
  expect_identical(class(out), class(sce))
  expect_true(methods::validObject(out))
  expect_identical(dimnames(out), dimnames(sce))
  expect_identical(
    SummarizedExperiment::assays(
      out,
      withDimnames = FALSE
    )[names(original_assays)],
    original_assays
  )
  expect_identical(
    SummarizedExperiment::rowData(out)[
      ,
      names(original_row_data),
      drop = FALSE
    ],
    original_row_data
  )
  expect_identical(
    SummarizedExperiment::colData(out)[
      ,
      names(original_col_data),
      drop = FALSE
    ],
    original_col_data
  )
  expect_identical(
    S4Vectors::metadata(out)[names(original_metadata)],
    original_metadata
  )
  expect_identical(
    SingleCellExperiment::reducedDim(
      out,
      "existing",
      withDimnames = FALSE
    ),
    original_reduced
  )
  expect_identical(
    SingleCellExperiment::altExp(out, "spike"),
    original_alt
  )
  expect_identical(
    SingleCellExperiment::rowPair(out, "existing_rows"),
    original_row_pair
  )
  expect_identical(
    SingleCellExperiment::colPair(out, "existing_cells"),
    original_col_pair
  )
  expect_identical(
    SingleCellExperiment::sizeFactors(out),
    original_size_factors
  )
  expect_identical(
    SingleCellExperiment::colLabels(out),
    original_labels
  )
  expect_identical(
    SingleCellExperiment::mainExpName(out),
    "endogenous"
  )

  expect_equal(
    SummarizedExperiment::assay(out, "cudacell_logcounts"),
    fit$normalized
  )
  expect_equal(
    SingleCellExperiment::reducedDim(out, "CUDACELL_PCA"),
    fit$pca$x
  )
  expect_identical(
    rownames(
      SingleCellExperiment::reducedDim(out, "CUDACELL_PCA")
    ),
    colnames(sce)
  )

  variable <- fit$variable_features
  row_data <- SummarizedExperiment::rowData(out)
  expect_equal(
    row_data$cudacell_mean[variable$index],
    variable$mean
  )
  expect_equal(
    row_data$cudacell_variance[variable$index],
    variable$variance
  )
  expect_equal(
    row_data$cudacell_dispersion[variable$index],
    variable$dispersion
  )
  expect_identical(
    row_data$cudacell_hvg[variable$index],
    variable$selected
  )
  expect_identical(
    row_data$cudacell_hvg_rank[variable$index],
    seq_len(nrow(sce))
  )

  hits <- as.data.frame(
    SingleCellExperiment::colPair(out, "CUDACELL_KNN")
  )
  hits <- hits[order(hits$from, hits$rank), , drop = FALSE]
  rownames(hits) <- NULL
  expect_identical(
    hits$from,
    rep(seq_len(ncol(sce)), each = arguments$k)
  )
  expect_identical(
    hits$to,
    as.integer(as.vector(t(fit$neighbors$index)))
  )
  expect_equal(
    hits$distance,
    as.numeric(as.vector(t(fit$neighbors$distance)))
  )
  expect_identical(
    hits$rank,
    rep(seq_len(arguments$k), times = ncol(sce))
  )

  expected_size_factors <- as.numeric(
    Matrix::colSums(SummarizedExperiment::assay(sce, "counts"))
  ) / 10000
  expect_equal(
    SummarizedExperiment::colData(out)$cudacell_size_factor,
    expected_size_factors
  )
  expect_identical(cuda_provenance(out), cuda_provenance(fit))
  expect_identical(
    S4Vectors::metadata(out)$cudacellr$schema,
    "cudacellr-sce/1"
  )
})

test_that("cudacell_sce maps duplicate feature names by position", {
  sce <- .sce_example()
  out <- cudacell_sce(
    sce,
    n_hvg = 5,
    n_components = 3,
    k = 3,
    device = "cpu"
  )
  fit <- cudacell_workflow(
    SummarizedExperiment::assay(sce, "counts"),
    n_hvg = 5,
    n_components = 3,
    k = 3,
    device = "cpu"
  )
  variable <- fit$variable_features

  expect_identical(rownames(out), rownames(sce))
  expect_identical(rownames(out)[1:2], rep("duplicated_gene", 2))
  expect_equal(
    SummarizedExperiment::rowData(out)$cudacell_dispersion[
      variable$index
    ],
    variable$dispersion
  )
})

test_that("SCE kNN colPairs remap automatically after cell subsetting", {
  sce <- .sce_example()
  out <- cudacell_sce(
    sce,
    n_hvg = 5,
    n_components = 3,
    k = 3,
    device = "cpu"
  )
  retained <- c(12L, 7L, 2L, 1L)
  original_hits <- as.data.frame(
    SingleCellExperiment::colPair(out, "CUDACELL_KNN")
  )
  expected <- original_hits[
    original_hits$from %in% retained &
      original_hits$to %in% retained,
    ,
    drop = FALSE
  ]
  expected$from <- match(expected$from, retained)
  expected$to <- match(expected$to, retained)
  expected <- expected[
    order(expected$from, expected$rank, expected$to),
    ,
    drop = FALSE
  ]
  rownames(expected) <- NULL

  subset_hits <- as.data.frame(
    SingleCellExperiment::colPair(
      out[, retained],
      "CUDACELL_KNN"
    )
  )
  subset_hits <- subset_hits[
    order(subset_hits$from, subset_hits$rank, subset_hits$to),
    ,
    drop = FALSE
  ]
  rownames(subset_hits) <- NULL
  expect_identical(subset_hits, expected)
})

test_that("output collisions fail before compute and overwrite is targeted", {
  sce <- .sce_example()
  out <- cudacell_sce(
    sce,
    n_hvg = 5,
    n_components = 3,
    k = 3,
    device = "cpu"
  )
  computed <- FALSE
  testthat::local_mocked_bindings(
    cudacell_workflow = function(...) {
      computed <<- TRUE
      stop("unexpected compute")
    },
    .package = "cudacellr"
  )
  condition <- tryCatch(
    cudacell_sce(
      out,
      n_hvg = 5,
      n_components = 3,
      k = 3,
      device = "cpu"
    ),
    error = identity
  )
  expect_s3_class(condition, "cudacell_sce_collision_error")
  expect_false(computed)
  expect_true(length(condition$conflicts) >= 5L)
  expect_match(condition$message, "overwrite = TRUE", fixed = TRUE)
})

test_that("overwrite replaces only named cudacell outputs", {
  sce <- .sce_example()
  out <- cudacell_sce(
    sce,
    n_hvg = 5,
    n_components = 3,
    k = 3,
    device = "cpu"
  )
  untouched <- SummarizedExperiment::assay(out, "untouched")
  existing_embedding <- SingleCellExperiment::reducedDim(
    out,
    "existing",
    withDimnames = FALSE
  )
  existing_metadata <- S4Vectors::metadata(out)$owner

  replaced <- cudacell_sce(
    out,
    n_hvg = 4,
    n_components = 2,
    k = 2,
    device = "cpu",
    overwrite = TRUE
  )
  expect_identical(
    SummarizedExperiment::assay(replaced, "untouched"),
    untouched
  )
  expect_identical(
    SingleCellExperiment::reducedDim(
      replaced,
      "existing",
      withDimnames = FALSE
    ),
    existing_embedding
  )
  expect_identical(
    S4Vectors::metadata(replaced)$owner,
    existing_metadata
  )
  expect_identical(
    ncol(
      SingleCellExperiment::reducedDim(
        replaced,
        "CUDACELL_PCA"
      )
    ),
    2L
  )
})

test_that("normalization controls and canonical size factors are explicit", {
  sce <- .sce_example()
  expect_error(
    cudacell_sce(
      sce,
      n_hvg = 5,
      n_components = 3,
      k = 3,
      device = "cpu",
      set_size_factors = TRUE
    ),
    class = "cudacell_sce_collision_error"
  )

  out <- cudacell_sce(
    sce,
    n_hvg = 5,
    n_components = 3,
    k = 3,
    device = "cpu",
    scale_factor = 1000,
    log1p = FALSE,
    scale. = FALSE,
    set_size_factors = TRUE,
    overwrite = TRUE
  )
  normalized <- SummarizedExperiment::assay(
    out,
    "cudacell_normalized"
  )
  expect_equal(
    as.numeric(Matrix::colSums(normalized)),
    rep(1000, ncol(out))
  )
  expect_identical(
    SingleCellExperiment::sizeFactors(out),
    SummarizedExperiment::colData(out)$cudacell_size_factor
  )
  expect_false(
    S4Vectors::metadata(out)$cudacellr$parameters$log1p
  )
  expect_false(
    S4Vectors::metadata(out)$cudacellr$parameters$scale
  )
})

test_that("adapter errors identify invalid SCE inputs and output names", {
  .skip_sce_dependencies()
  expect_error(
    cudacell_sce(matrix(1, 3, 3)),
    class = "cudacell_sce_type_error"
  )
  sce <- .sce_example()
  expect_error(
    cudacell_sce(sce, assay = "missing"),
    class = "cudacell_sce_assay_error"
  )
  expect_error(
    cudacell_sce(sce, normalized_assay = ""),
    class = "cudacell_sce_name_error"
  )
  expect_error(
    cudacell_sce(
      sce,
      normalized_assay = "counts",
      overwrite = TRUE
    ),
    class = "cudacell_sce_name_error"
  )
  expect_error(
    cudacell_sce(sce, overwrite = NA),
    class = "cudacell_sce_name_error"
  )

  invalid <- sce
  SummarizedExperiment::assay(
    invalid,
    "character_counts"
  ) <- matrix(
    "x",
    nrow = nrow(sce),
    ncol = ncol(sce),
    dimnames = dimnames(sce)
  )
  expect_error(
    cudacell_sce(invalid, assay = "character_counts"),
    "numeric values"
  )
})

test_that("delayed assays require explicit realization", {
  .skip_sce_dependencies()
  skip_if_not_installed("DelayedArray")
  sce <- .sce_example()
  delayed <- DelayedArray::DelayedArray(
    as.matrix(SummarizedExperiment::assay(sce, "counts"))
  )
  delayed_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = delayed)
  )
  expect_error(
    cudacell_sce(
      delayed_sce,
      n_hvg = 5,
      n_components = 3,
      k = 3,
      device = "cpu"
    ),
    class = "cudacell_sce_assay_error"
  )

  out <- cudacell_sce(
    delayed_sce,
    n_hvg = 5,
    n_components = 3,
    k = 3,
    device = "cpu",
    realize = TRUE
  )
  expect_s4_class(out, "SingleCellExperiment")
  expect_true(
    S4Vectors::metadata(out)$cudacellr$input_assay_materialized
  )
  expect_true(
    "DelayedMatrix" %in%
      S4Vectors::metadata(out)$cudacellr$input_assay_class
  )
})

test_that("strict CUDA validation precedes delayed assay realization", {
  .skip_sce_dependencies()
  skip_if_not_installed("DelayedArray")
  counts <- DelayedArray::DelayedArray(
    matrix(1, nrow = 6, ncol = 6)
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )
  unavailable <- structure(
    list(
      torch_installed = FALSE,
      torch_version = NA_character_,
      cuda_available = FALSE,
      cuda_device_count = 0L,
      reason = "torch_not_installed",
      detection_error = NULL
    ),
    class = "cuda_diagnostics"
  )
  testthat::local_mocked_bindings(
    cuda_diagnostics = function() unavailable,
    .package = "cudaverse"
  )

  condition <- tryCatch(
    cudacell_sce(
      sce,
      n_hvg = 5,
      n_components = 3,
      k = 3,
      device = "cuda"
    ),
    error = identity
  )
  expect_s3_class(condition, "cudaverse_cuda_unavailable")
})

test_that("objects without dimnames retain NULL identifiers", {
  sce <- .sce_example(with_names = FALSE)
  out <- cudacell_sce(
    sce,
    n_hvg = 5,
    n_components = 3,
    k = 3,
    device = "cpu"
  )

  expect_null(rownames(out))
  expect_null(colnames(out))
  expect_null(
    rownames(
      SummarizedExperiment::assay(out, "cudacell_logcounts")
    )
  )
  expect_null(
    rownames(
      SingleCellExperiment::reducedDim(out, "CUDACELL_PCA")
    )
  )
})

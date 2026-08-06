.cell_test_counts <- function() {
  counts <- matrix(
    c(
      4, 0, 1, 3, 2, 0,
      0, 2, 3, 1, 0, 4,
      5, 1, 0, 2, 3, 1,
      1, 4, 2, 0, 1, 3,
      2, 1, 4, 3, 0, 2,
      3, 2, 1, 4, 2, 0
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(
      paste0("gene_", seq_len(6)),
      paste0("cell_", seq_len(6))
    )
  )
  Matrix::Matrix(counts, sparse = TRUE)
}

test_that("single-cell stages compose into one inspectable workflow", {
  counts <- .cell_test_counts()
  normalized <- cuda_normalize_counts(counts)
  expect_identical(
    cuda_provenance(normalized)$stage,
    "normalization"
  )

  variable <- cuda_hvg(normalized, n_top = 4)
  expect_identical(
    cuda_provenance(variable)$stage,
    c("normalization", "hvg")
  )

  pca <- cuda_cell_pca(
    counts,
    n_components = 2,
    n_hvg = 4,
    device = "cpu"
  )
  expect_identical(
    cuda_provenance(pca)$stage,
    c(
      "normalization",
      "hvg",
      "pca_preprocessing",
      "pca_decomposition"
    )
  )
  expect_identical(pca$compute_device, "cpu")

  neighbors <- cuda_cell_neighbors(
    pca,
    k = 2,
    batch_size = 3,
    device = "cpu"
  )
  expect_identical(
    cuda_provenance(neighbors)$stage,
    c(
      "normalization",
      "hvg",
      "pca_preprocessing",
      "pca_decomposition",
      "knn_distance",
      "knn_neighbor_selection"
    )
  )

  workflow <- cudacell_workflow(
    counts,
    n_hvg = 4,
    n_components = 2,
    k = 2,
    batch_size = 3,
    device = "cpu"
  )
  expect_identical(
    cuda_provenance(workflow)$stage,
    cuda_provenance(neighbors)$stage
  )
  expect_identical(workflow$compute_device, "cpu")
  expect_output(print(workflow), "compute=cpu")
})

test_that("upstream sparse provenance is retained", {
  sparse <- cudaverse::cuda_sparse(
    .cell_test_counts(),
    device = "cpu"
  )
  fit <- cuda_cell_pca(
    sparse,
    n_components = 2,
    n_hvg = 4,
    device = "cpu"
  )
  provenance <- cuda_provenance(fit)

  expect_identical(
    provenance$stage[[1L]],
    "source_sparse_materialization"
  )
  expect_true(all(
    c("normalization", "hvg", "pca_decomposition") %in%
      provenance$stage
  ))
})

test_that("declared but invalid upstream provenance is not discarded", {
  counts <- .cell_test_counts()
  attr(counts, "compute_stages") <- list(
    source = cudaverse::cuda_stage(
      requested_device = "fixed-cpu",
      device = "cpu",
      backend = "Matrix",
      selection_reason = "algorithm_cpu_only"
    )
  )
  attr(counts, "compute_device") <- "cpu"
  attr(counts, "provenance_schema") <- "cudaverse-stage/99"

  condition <- tryCatch(
    cuda_normalize_counts(counts),
    error = identity
  )
  expect_s3_class(condition, "cudaverse_provenance_schema_error")

  attr(counts, "provenance_schema") <- "cudaverse-stage/1"
  attr(counts, "compute_stages") <- list(
    broken = list(device = "cpu")
  )
  expect_error(
    cuda_normalize_counts(counts),
    "does not follow the cudaverse stage schema"
  )
})

test_that("workflow validates strict CUDA before CPU preprocessing", {
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
    cudacell_workflow("not a matrix", device = "cuda"),
    error = identity
  )
  expect_s3_class(condition, "cudaverse_cuda_unavailable")

  fit <- cudacell_workflow(
    .cell_test_counts(),
    n_hvg = 4,
    n_components = 2,
    k = 2,
    device = "auto"
  )
  provenance <- cuda_provenance(fit)
  requested <- provenance$requested_device == "auto"
  expect_true(any(requested))
  expect_true(all(provenance$fallback[requested]))
  expect_true(all(
    provenance$selection_reason[requested] == "torch_not_installed"
  ))
})

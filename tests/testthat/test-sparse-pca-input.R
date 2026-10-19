.sparse_pca_counts <- function() {
  set.seed(1)
  x <- matrix(rpois(30 * 20, lambda = 2), 30, 20,
              dimnames = list(paste0("gene", 1:30), paste0("cell", 1:20)))
  Matrix::Matrix(x, sparse = TRUE)
}

test_that("native preparation preserves sparse feature positions and names", {
  counts <- .sparse_pca_counts()
  rownames(counts)[1:2] <- "duplicate"
  normalized <- cuda_normalize_counts(counts)
  indices <- c(2L, 5L, 1L, 8L)
  features <- rownames(counts)[indices]
  cpu <- .cell_pca_input(normalized, indices, features, "cpu")
  select <- cudaverse::cuda_select_device
  testthat::local_mocked_bindings(
    cuda_select_device = function(device = c("auto", "cuda", "cpu")) {
      if (identical(device, "cpu")) return(select("cpu"))
      list(device = "cuda", backend = "native")
    },
    .package = "cudaverse"
  )
  for (device in c("cuda", "auto")) {
    input <- .cell_pca_input(normalized, indices, features, device)
    expect_s3_class(input, "cudasparse")
    expect_identical(input$device, "cpu")
    expect_equal(as.matrix(cudaverse::to_dgCMatrix(input)), cpu)
    expect_identical(input$dimnames, dimnames(cpu))
  }
})

test_that("torch preparation retains the dense contract", {
  normalized <- cuda_normalize_counts(.sparse_pca_counts())
  testthat::local_mocked_bindings(
    cuda_select_device = function(...) list(device = "cuda", backend = "torch"),
    .package = "cudaverse"
  )
  input <- .cell_pca_input(normalized, 1:4, rownames(normalized)[1:4], "cuda")
  expect_true(is.matrix(input))
  expect_equal(input, t(as.matrix(normalized[1:4, ])))
})

test_that("invalid component requests fail before input preparation", {
  counts <- .sparse_pca_counts()
  variable <- cuda_hvg(counts, n_top = 4)
  testthat::local_mocked_bindings(
    .cell_pca_input = function(...) stop("Input was allocated"),
    .package = "cudacellr"
  )
  expect_error(.cell_pca_from_normalized(counts, variable, 5, TRUE, "cpu"),
               "cannot exceed 4")
  expect_error(.cell_pca_from_normalized(counts, variable, 2^31, TRUE, "cpu"),
               "positive whole number")
})

test_that("native sparse preparation preserves PCA, ties and score reuse", {
  selection <- cudaverse::cuda_select_device("auto")
  available <- identical(selection$device, "cuda") &&
    identical(selection$backend, "native")
  if (identical(Sys.getenv("CUDACELLR_REQUIRE_NATIVE"), "true")) {
    expect_true(available)
    if (!available) stop("Native CUDA hardware validation was required")
  } else {
    skip_if_not(available, "Native CUDA hardware unavailable")
  }
  counts <- .sparse_pca_counts()
  counts[, 2] <- counts[, 1]
  rownames(counts)[1:2] <- "duplicate"
  normalized <- cuda_normalize_counts(counts)
  variable <- cuda_hvg(normalized, n_top = 10)
  selected <- variable[variable$selected, , drop = FALSE]
  for (scaled in c(FALSE, TRUE)) {
    reference <- cudaverse::cuda_pca(
      t(as.matrix(normalized[selected$index, , drop = FALSE])),
      n_components = 3, scale. = scaled, device = "cuda"
    )
    candidate <- .cell_pca_from_normalized(
      normalized, variable, 3, scaled, "cuda"
    )
    expect_identical(dimnames(candidate$x), dimnames(reference$x))
    expect_identical(rownames(candidate$rotation), selected$feature)
    expect_equal(candidate$sdev, reference$sdev, tolerance = 1e-10)
    expect_equal(tcrossprod(candidate$rotation), tcrossprod(reference$rotation),
                 tolerance = 1e-8)
    expect_equal(tcrossprod(candidate$x, candidate$rotation),
                 tcrossprod(reference$x, reference$rotation), tolerance = 1e-8)
    expect_identical(as.numeric(candidate$x[1, ]), as.numeric(candidate$x[2, ]))
    a <- cuda_cell_neighbors(reference, k = 4, device = "cuda", batch_size = 3)
    b <- cuda_cell_neighbors(candidate, k = 4, device = "cuda", batch_size = 7)
    expect_identical(b$index, a$index)
    expect_equal(b$distance, a$distance, tolerance = 1e-8)
    stages <- cuda_provenance(b)
    expect_true(any(stages$backend == "native"))
    input_stage <- stages[stages$stage == "knn_input_materialization", ]
    expect_identical(input_stage$selection_reason, "device_resident_input")
    expect_identical(input_stage$output_device, "cuda")
  }
})

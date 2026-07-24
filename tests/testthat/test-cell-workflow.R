example_counts <- function() {
  set.seed(1)
  counts <- matrix(rpois(30 * 20, lambda = 2), 30, 20)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  Matrix::Matrix(counts, sparse = TRUE)
}

test_that("normalization preserves sparse shape and library scaling", {
  normalized <- cuda_normalize_counts(
    example_counts(), scale_factor = 1000, log1p = FALSE
  )

  expect_s4_class(normalized, "dgCMatrix")
  expect_identical(dim(normalized), c(30L, 20L))
  expect_equal(as.numeric(Matrix::colSums(normalized)), rep(1000, 20))
})

test_that("HVG selection returns ranked feature statistics", {
  hvg <- cuda_hvg(example_counts(), n_top = 8)

  expect_equal(sum(hvg$selected), 8)
  expect_true(all(diff(hvg$dispersion) <= 0))
  expect_true(all(c("mean", "variance", "dispersion") %in% names(hvg)))
})

test_that("PCA and neighbors compose", {
  pca <- cuda_cell_pca(
    example_counts(), n_components = 3, n_hvg = 10, device = "cpu"
  )
  neighbours <- cuda_cell_neighbors(pca, k = 4, device = "cpu")

  expect_s3_class(pca, "cuda_pca")
  expect_identical(dim(pca$x), c(20L, 3L))
  expect_length(pca$features, 10)
  expect_identical(dim(neighbours$index), c(20L, 4L))
})

test_that("workflow returns all analysis stages", {
  fit <- cudacell_workflow(
    example_counts(),
    n_hvg = 10,
    n_components = 3,
    k = 4,
    device = "cpu"
  )
  expect_s3_class(fit, "cudacell_workflow")
  expect_named(
    fit,
    c("normalized", "variable_features", "pca", "neighbors")
  )
})

test_that("invalid count matrices fail clearly", {
  counts <- as.matrix(example_counts())
  counts[[1]] <- -1
  expect_error(cuda_normalize_counts(counts), "non-negative")

  counts <- as.matrix(example_counts())
  counts[, 1] <- 0
  expect_error(cuda_normalize_counts(counts), "positive library")

  expect_error(cuda_cell_pca(counts, n_hvg = NA), "positive whole")
  expect_error(cuda_cell_pca(counts, n_components = 0), "positive whole")
  expect_error(cudacell_workflow(counts, n_hvg = 1.5), "positive whole")
})

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
  expect_identical(dimnames(normalized), dimnames(example_counts()))
  expect_equal(as.numeric(Matrix::colSums(normalized)), rep(1000, 20))
})

test_that("HVG selection returns ranked feature statistics", {
  hvg <- cuda_hvg(example_counts(), n_top = 8)

  expect_equal(sum(hvg$selected), 8)
  expect_true(all(diff(hvg$dispersion) <= 0))
  expect_true(all(
    c("index", "feature", "mean", "variance", "dispersion") %in% names(hvg)
  ))
})

test_that("PCA and neighbors compose", {
  pca <- cuda_cell_pca(
    example_counts(), n_components = 3, n_hvg = 10, device = "cpu"
  )
  neighbours <- cuda_cell_neighbors(
    pca,
    k = 4,
    device = "cpu",
    batch_size = 2
  )
  one_at_a_time <- cuda_cell_neighbors(
    pca,
    k = 4,
    device = "cpu",
    batch_size = 1
  )

  expect_s3_class(pca, "cuda_pca")
  expect_identical(dim(pca$x), c(20L, 3L))
  expect_length(pca$features, 10)
  expect_identical(rownames(pca$x), colnames(example_counts()))
  expect_identical(rownames(pca$rotation), pca$features)
  expect_identical(rownames(neighbours$index), colnames(example_counts()))
  expect_identical(dimnames(neighbours$distance), dimnames(neighbours$index))
  expect_identical(dim(neighbours$index), c(20L, 4L))
  expect_identical(neighbours$index, one_at_a_time$index)
  expect_equal(neighbours$distance, one_at_a_time$distance)
})

test_that("workflow returns all analysis stages", {
  fit <- cudacell_workflow(
    example_counts(),
    n_hvg = 10,
    n_components = 3,
    k = 4,
    batch_size = 2,
    device = "cpu"
  )
  expect_s3_class(fit, "cudacell_workflow")
  expect_named(
    fit,
    c("normalized", "variable_features", "pca", "neighbors")
  )
  expect_identical(rownames(fit$pca$x), colnames(example_counts()))
  expect_identical(
    rownames(fit$neighbors$index),
    colnames(example_counts())
  )
  expect_output(print(fit), "<cudacell_workflow")
})

test_that("workflow reuses normalization and HVG preprocessing", {
  normalizations <- 0L
  selections <- 0L
  normalize_original <- cuda_normalize_counts
  hvg_original <- cuda_hvg
  testthat::local_mocked_bindings(
    cuda_normalize_counts = function(...) {
      normalizations <<- normalizations + 1L
      normalize_original(...)
    },
    cuda_hvg = function(...) {
      selections <<- selections + 1L
      hvg_original(...)
    },
    .package = "cudacellr"
  )

  fit <- cudacell_workflow(
    example_counts(),
    n_hvg = 10,
    n_components = 3,
    k = 4,
    device = "cpu"
  )

  expect_s3_class(fit, "cudacell_workflow")
  expect_identical(normalizations, 1L)
  expect_identical(selections, 1L)
  expect_setequal(
    fit$pca$features,
    fit$variable_features$feature[fit$variable_features$selected]
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
  expect_error(cuda_normalize_counts(example_counts(), log1p = NA), "TRUE or FALSE")
})

test_that("duplicate feature names retain their original positions", {
  counts <- example_counts()
  rownames(counts)[1:2] <- "duplicated_gene"
  variable <- cuda_hvg(counts, n_top = 10)
  pca <- cudacellr:::.cell_pca_from_normalized(
    normalized = cuda_normalize_counts(counts),
    variable = variable,
    n_components = 3,
    scale. = TRUE,
    device = "cpu"
  )

  selected <- variable[variable$selected, , drop = FALSE]
  expect_identical(pca$features, selected$feature)
  expect_identical(rownames(pca$rotation), selected$feature)
})

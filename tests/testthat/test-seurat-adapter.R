.skip_seurat_dependency <- function() {
  skip_if_not_installed("SeuratObject", minimum_version = "5.0.0")
}

.seurat_fixture_tool <- function(object) {
  SeuratObject::Tool(object) <- list(
    owner = "preserve-me",
    nested = list(value = 42L)
  )
  object
}

.seurat_example <- function() {
  .skip_seurat_dependency()
  old <- options(Seurat.object.assay.version = "v5")
  on.exit(options(old), add = TRUE)

  counts <- matrix(
    c(
      4, 0, 1, 3, 2, 1, 0, 2, 1, 4, 2, 3,
      0, 2, 3, 1, 1, 4, 2, 0, 3, 1, 2, 4,
      5, 1, 0, 2, 3, 1, 4, 2, 1, 0, 3, 2,
      1, 4, 2, 1, 1, 3, 2, 4, 0, 2, 1, 3,
      2, 1, 4, 3, 1, 2, 0, 1, 4, 3, 2, 1,
      3, 2, 1, 4, 2, 1, 3, 0, 2, 4, 1, 2,
      1, 3, 2, 4, 0, 2, 1, 3, 2, 1, 4, 2,
      2, 4, 1, 0, 3, 2, 4, 1, 2, 3, 1, 4
    ),
    nrow = 8,
    byrow = TRUE,
    dimnames = list(
      paste0("gene", seq_len(8)),
      paste0("cell_", seq_len(12))
    )
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  cell_metadata <- data.frame(
    batch = factor(rep(c("one", "two"), each = 6)),
    score = seq_len(ncol(counts)),
    row.names = colnames(counts)
  )
  object <- SeuratObject::CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    meta.data = cell_metadata,
    project = "preserve-project"
  )

  second <- counts[seq_len(3), , drop = FALSE]
  second@x <- second@x + 10
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(
    counts = second
  )

  existing_embedding <- matrix(
    seq_len(ncol(counts) * 2L),
    nrow = ncol(counts),
    dimnames = list(
      colnames(counts),
      c("EXISTING_1", "EXISTING_2")
    )
  )
  object[["existing_reduction"]] <-
    SeuratObject::CreateDimReducObject(
      embeddings = existing_embedding,
      assay = "RNA",
      key = "EXISTING_"
    )

  existing_index <- cbind(
    c(2:12, 1),
    c(3:12, 1, 2)
  )
  rownames(existing_index) <- colnames(counts)
  existing_distance <- matrix(
    seq_len(length(existing_index)) / 100,
    nrow = ncol(counts),
    dimnames = dimnames(existing_index)
  )
  object[["existing_neighbor"]] <- methods::new(
    "Neighbor",
    nn.idx = existing_index,
    nn.dist = existing_distance,
    alg.idx = NULL,
    alg.info = list(source = "fixture"),
    cell.names = colnames(counts)
  )

  adjacency <- Matrix::sparseMatrix(
    i = c(seq_len(11), 2:12),
    j = c(2:12, seq_len(11)),
    x = 1,
    dims = rep(ncol(counts), 2),
    dimnames = list(colnames(counts), colnames(counts))
  )
  object[["existing_graph"]] <- SeuratObject::as.Graph(adjacency)
  SeuratObject::DefaultAssay(object) <- "RNA"
  SeuratObject::Idents(object) <- factor(
    rep(c("A", "B"), each = 6)
  )
  SeuratObject::Misc(object, slot = "owner") <- list(
    value = "preserve-me"
  )
  .seurat_fixture_tool(object)
}

.seurat_workflow_args <- function() {
  list(
    n_hvg = 6,
    n_components = 3,
    k = 3,
    batch_size = 2,
    device = "cpu"
  )
}

test_that("cudacell_seurat writes native outputs and preserves object state", {
  object <- .seurat_example()
  original <- object
  original_assays <- lapply(
    SeuratObject::Assays(object),
    function(name) object[[name]]
  )
  names(original_assays) <- SeuratObject::Assays(object)
  original_reduction <- object[["existing_reduction"]]
  original_neighbor <- object[["existing_neighbor"]]
  original_graph <- object[["existing_graph"]]
  original_metadata <- object[[]]
  original_idents <- SeuratObject::Idents(object)
  original_misc <- SeuratObject::Misc(object)
  original_tools <- SeuratObject::Tool(object)

  arguments <- .seurat_workflow_args()
  out <- do.call(
    cudacell_seurat,
    c(list(x = object), arguments)
  )
  counts <- SeuratObject::LayerData(
    object[["RNA"]],
    layer = "counts",
    fast = FALSE
  )
  fit <- do.call(
    cudacell_workflow,
    c(list(counts = counts), arguments)
  )

  expect_identical(object, original)
  expect_s4_class(out, "Seurat")
  expect_true(methods::validObject(out))
  expect_identical(
    SeuratObject::DefaultAssay(out),
    SeuratObject::DefaultAssay(object)
  )
  expect_identical(SeuratObject::Idents(out), original_idents)
  expect_identical(
    SeuratObject::Project(out),
    "preserve-project"
  )
  for (name in names(original_assays)) {
    expect_identical(out[[name]], original_assays[[name]])
  }
  expect_identical(out[["existing_reduction"]], original_reduction)
  expect_identical(out[["existing_neighbor"]], original_neighbor)
  expect_identical(out[["existing_graph"]], original_graph)
  expect_identical(
    out[[]][, names(original_metadata), drop = FALSE],
    original_metadata
  )
  expect_identical(
    SeuratObject::Misc(out),
    original_misc
  )
  expect_identical(
    SeuratObject::Tool(out)[
      SeuratObject::Tool(out) %in% original_tools
    ],
    original_tools
  )

  expect_true("CUDACELL" %in% SeuratObject::Assays(out))
  expect_s4_class(out[["CUDACELL"]], "Assay5")
  expect_identical(
    SeuratObject::Layers(out[["CUDACELL"]]),
    "data"
  )
  normalized <- SeuratObject::LayerData(
    out[["CUDACELL"]],
    layer = "data",
    fast = FALSE
  )
  expect_equal(as.matrix(normalized), as.matrix(fit$normalized))
  expect_identical(dimnames(normalized), dimnames(counts))

  feature_metadata <- out[["CUDACELL"]][[]]
  variable <- fit$variable_features
  expect_equal(
    feature_metadata$cudacell_mean[variable$index],
    variable$mean
  )
  expect_equal(
    feature_metadata$cudacell_variance[variable$index],
    variable$variance
  )
  expect_equal(
    feature_metadata$cudacell_dispersion[variable$index],
    variable$dispersion
  )
  expect_identical(
    feature_metadata$cudacell_hvg[variable$index],
    variable$selected
  )
  expect_identical(
    feature_metadata$cudacell_hvg_rank[variable$index],
    seq_len(nrow(counts))
  )
  expect_identical(
    order(feature_metadata$cudacell_hvg_rank),
    variable$index
  )
  expect_identical(
    feature_metadata$cudacell_hvg,
    feature_metadata$cudacell_hvg_rank <= arguments$n_hvg
  )
  expect_identical(
    SeuratObject::VariableFeatures(out[["CUDACELL"]]),
    rownames(counts)[variable$index[variable$selected]]
  )

  reduction <- out[["cudacell_pca"]]
  expect_s4_class(reduction, "DimReduc")
  expect_equal(
    unname(SeuratObject::Embeddings(reduction)),
    unname(fit$pca$x)
  )
  expect_identical(
    rownames(SeuratObject::Embeddings(reduction)),
    colnames(counts)
  )
  expect_equal(
    unname(SeuratObject::Loadings(reduction)),
    unname(fit$pca$rotation)
  )
  expect_equal(SeuratObject::Stdev(reduction), fit$pca$sdev)
  expect_identical(
    SeuratObject::DefaultAssay(reduction),
    "CUDACELL"
  )
  expect_identical(
    SeuratObject::Key(reduction),
    "CUDACELLPC_"
  )

  neighbor <- out[["cudacell_knn"]]
  expect_s4_class(neighbor, "Neighbor")
  expect_identical(
    SeuratObject::Indices(neighbor),
    fit$neighbors$index
  )
  expect_equal(
    SeuratObject::Distances(neighbor),
    fit$neighbors$distance
  )
  expect_identical(
    SeuratObject::Cells(neighbor),
    colnames(counts)
  )

  expected_size_factors <- as.numeric(Matrix::colSums(counts)) / 10000
  expect_equal(
    unname(out$cudacell_size_factor),
    expected_size_factors
  )
  tool <- SeuratObject::Tool(
    out,
    slot = "cudacell_seurat"
  )
  expect_identical(tool$schema, "cudacellr-seurat/1")
  expect_identical(tool$input_assay, "RNA")
  expect_identical(tool$input_layer, "counts")
  expect_false(tool$input_layer_materialized)
  expect_identical(tool$outputs$assay, "CUDACELL")
  expect_identical(tool$outputs$reduction, "cudacell_pca")
  expect_identical(tool$outputs$neighbor, "cudacell_knn")
  expect_identical(cuda_provenance(out), cuda_provenance(fit))
  expect_identical(
    cudaverse::cuda_provenance(out),
    cuda_provenance(fit)
  )
})

test_that("custom output names remain native and fully namespaced", {
  object <- .seurat_example()
  out <- cudacell_seurat(
    object,
    n_hvg = 5,
    n_components = 2,
    k = 2,
    device = "cpu",
    output_assay = "CUDARESULT",
    reduction = "cuda_result_pca",
    reduction_key = "CUDARESULTPC_",
    neighbor = "cuda_result_knn",
    metadata_prefix = "cudaresult"
  )

  expect_s4_class(out[["CUDARESULT"]], "Assay5")
  expect_s4_class(out[["cuda_result_pca"]], "DimReduc")
  expect_s4_class(out[["cuda_result_knn"]], "Neighbor")
  expect_true("cudaresult_size_factor" %in% colnames(out[[]]))
  expect_true(
    all(
      paste0(
        "cudaresult_",
        c("mean", "variance", "dispersion", "hvg", "hvg_rank")
      ) %in% colnames(out[["CUDARESULT"]][[]])
    )
  )
  tool <- SeuratObject::Tool(
    out,
    slot = "cudacell_seurat"
  )
  expect_identical(tool$outputs$assay, "CUDARESULT")
  expect_identical(tool$outputs$reduction_key, "CUDARESULTPC_")
})

test_that("output collisions fail before compute", {
  object <- .seurat_example()
  out <- cudacell_seurat(
    object,
    n_hvg = 6,
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
    cudacell_seurat(
      out,
      n_hvg = 6,
      n_components = 3,
      k = 3,
      device = "cpu"
    ),
    error = identity
  )
  expect_s3_class(
    condition,
    "cudacell_seurat_collision_error"
  )
  expect_false(computed)
  expect_true(length(condition$conflicts) >= 5L)
  expect_match(condition$message, "overwrite = TRUE", fixed = TRUE)
})

test_that("overwrite replaces only named cudacell outputs", {
  object <- .seurat_example()
  first <- cudacell_seurat(
    object,
    n_hvg = 6,
    n_components = 3,
    k = 3,
    device = "cpu"
  )
  first_copy <- first
  original_rna <- first[["RNA"]]
  original_adt <- first[["ADT"]]
  original_reduction <- first[["existing_reduction"]]
  original_neighbor <- first[["existing_neighbor"]]
  original_graph <- first[["existing_graph"]]
  original_idents <- SeuratObject::Idents(first)
  original_misc <- SeuratObject::Misc(first)
  original_fixture_tool <- SeuratObject::Tool(
    first,
    slot = ".seurat_fixture_tool"
  )

  replaced <- cudacell_seurat(
    first,
    n_hvg = 5,
    n_components = 2,
    k = 2,
    device = "cpu",
    scale_factor = 1000,
    log1p = FALSE,
    scale. = FALSE,
    overwrite = TRUE
  )

  expect_identical(first, first_copy)
  expect_identical(replaced[["RNA"]], original_rna)
  expect_identical(replaced[["ADT"]], original_adt)
  expect_identical(
    replaced[["existing_reduction"]],
    original_reduction
  )
  expect_identical(
    replaced[["existing_neighbor"]],
    original_neighbor
  )
  expect_identical(replaced[["existing_graph"]], original_graph)
  expect_identical(SeuratObject::Idents(replaced), original_idents)
  expect_identical(SeuratObject::Misc(replaced), original_misc)
  expect_identical(
    SeuratObject::Tool(
      replaced,
      slot = ".seurat_fixture_tool"
    ),
    original_fixture_tool
  )
  expect_identical(
    ncol(SeuratObject::Embeddings(replaced[["cudacell_pca"]])),
    2L
  )
  expect_identical(
    ncol(SeuratObject::Indices(replaced[["cudacell_knn"]])),
    2L
  )
  expect_equal(
    unname(Matrix::colSums(
      SeuratObject::LayerData(
        replaced[["CUDACELL"]],
        layer = "data",
        fast = FALSE
      )
    )),
    rep(1000, ncol(replaced))
  )
  tool <- SeuratObject::Tool(
    replaced,
    slot = "cudacell_seurat"
  )
  expect_false(tool$parameters$log1p)
  expect_false(tool$parameters$scale)
})

test_that("unrelated subobject and key collisions are never overwritten", {
  object <- .seurat_example()
  object[["CUDACELL"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      rnorm(ncol(object) * 2),
      nrow = ncol(object),
      dimnames = list(
        colnames(object),
        c("OTHER_1", "OTHER_2")
      )
    ),
    assay = "RNA",
    key = "OTHER_"
  )
  expect_error(
    cudacell_seurat(
      object,
      n_hvg = 6,
      n_components = 3,
      k = 3,
      device = "cpu",
      overwrite = TRUE
    ),
    class = "cudacell_seurat_collision_error"
  )

  expect_error(
    cudacell_seurat(
      .seurat_example(),
      n_hvg = 6,
      n_components = 3,
      k = 3,
      device = "cpu",
      reduction_key = "cudacell_",
      overwrite = TRUE
    ),
    class = "cudacell_seurat_collision_error"
  )

  object <- .seurat_example()
  object[["other_reduction"]] <-
    SeuratObject::CreateDimReducObject(
      embeddings = matrix(
        rnorm(ncol(object) * 2),
        nrow = ncol(object),
        dimnames = list(
          colnames(object),
          c("CUDACELLPC_1", "CUDACELLPC_2")
        )
      ),
      assay = "RNA",
      key = "CUDACELLPC_"
    )
  expect_error(
    cudacell_seurat(
      object,
      n_hvg = 6,
      n_components = 3,
      k = 3,
      device = "cpu",
      overwrite = TRUE
    ),
    class = "cudacell_seurat_collision_error"
  )

  object <- .seurat_example()
  other_assay <- SeuratObject::CreateAssay5Object(
    counts = SeuratObject::LayerData(
      object[["RNA"]],
      layer = "counts",
      fast = FALSE
    )
  )
  SeuratObject::Key(other_assay) <- "CUDACELLPC_"
  object[["other_key_assay"]] <- other_assay
  expect_error(
    cudacell_seurat(
      object,
      n_hvg = 6,
      n_components = 3,
      k = 3,
      device = "cpu",
      overwrite = TRUE
    ),
    class = "cudacell_seurat_collision_error"
  )

  object <- .seurat_example()
  object[["other_key_reduction"]] <-
    SeuratObject::CreateDimReducObject(
      embeddings = matrix(
        rnorm(ncol(object) * 2),
        nrow = ncol(object),
        dimnames = list(
          colnames(object),
          c("cudacell_1", "cudacell_2")
        )
      ),
      assay = "RNA",
      key = "cudacell_"
    )
  expect_error(
    cudacell_seurat(
      object,
      n_hvg = 6,
      n_components = 3,
      k = 3,
      device = "cpu",
      overwrite = TRUE
    ),
    class = "cudacell_seurat_collision_error"
  )
})

test_that("adapter validates exact assay, layer, and output names", {
  .skip_seurat_dependency()
  expect_error(
    cudacell_seurat(matrix(1, 3, 3)),
    class = "cudacell_seurat_type_error"
  )
  object <- .seurat_example()
  expect_error(
    cudacell_seurat(object, assay = "missing"),
    class = "cudacell_seurat_layer_error"
  )
  expect_error(
    cudacell_seurat(object, layer = "missing"),
    class = "cudacell_seurat_layer_error"
  )
  expect_error(
    cudacell_seurat(object, output_assay = "RNA"),
    class = "cudacell_seurat_name_error"
  )
  expect_error(
    cudacell_seurat(
      object,
      output_assay = "same",
      reduction = "same"
    ),
    class = "cudacell_seurat_name_error"
  )
  expect_error(
    cudacell_seurat(object, reduction_key = "invalid-key"),
    class = "cudacell_seurat_name_error"
  )
  expect_error(
    cudacell_seurat(object, overwrite = NA),
    class = "cudacell_seurat_name_error"
  )
  expect_error(
    cudacell_seurat(object, realize = 1),
    class = "cudacell_seurat_name_error"
  )
})

test_that("selected layer must cover every object cell in order", {
  object <- .seurat_example()
  partial <- SeuratObject::LayerData(
    object[["RNA"]],
    layer = "counts",
    fast = FALSE
  )[, seq_len(8), drop = FALSE]
  object[["PARTIAL"]] <- SeuratObject::CreateAssay5Object(
    counts = partial
  )
  expect_error(
    cudacell_seurat(
      object,
      assay = "PARTIAL",
      n_hvg = 5,
      n_components = 3,
      k = 3,
      device = "cpu"
    ),
    class = "cudacell_seurat_mapping_error"
  )
})

test_that("non-memory-backed layers require explicit realization", {
  values <- matrix(
    seq_len(24),
    nrow = 4,
    dimnames = list(
      paste0("gene_", 1:4),
      paste0("cell_", 1:6)
    )
  )
  fake <- as.data.frame(values)
  expect_error(
    cudacellr:::.seurat_realize_layer(
      fake,
      assay = "RNA",
      layer = "counts",
      realize = FALSE
    ),
    class = "cudacell_seurat_layer_error"
  )
  realized <- cudacellr:::.seurat_realize_layer(
    fake,
    assay = "RNA",
    layer = "counts",
    realize = TRUE
  )
  expect_true(realized$materialized)
  expect_s4_class(realized$counts, "Matrix")
  expect_equal(as.matrix(realized$counts), values)
  expect_identical(dimnames(realized$counts), dimnames(values))
  expect_identical(realized$class, "data.frame")
})

test_that("strict CUDA validation precedes reading or realizing a layer", {
  object <- .seurat_example()
  layer_read <- FALSE
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
    .seurat_layer = function(...) {
      layer_read <<- TRUE
      stop("unexpected layer read")
    },
    .package = "cudacellr"
  )
  testthat::local_mocked_bindings(
    cuda_diagnostics = function() unavailable,
    .package = "cudaverse"
  )

  condition <- tryCatch(
    cudacell_seurat(
      object,
      n_hvg = 6,
      n_components = 3,
      k = 3,
      device = "cuda",
      realize = TRUE
    ),
    error = identity
  )
  expect_s3_class(condition, "cudaverse_cuda_unavailable")
  expect_false(layer_read)
})

test_that("Seurat provenance rejects missing and corrupt tool records", {
  object <- .seurat_example()
  expect_error(
    cuda_provenance(object),
    class = "cudacell_seurat_metadata_error"
  )

  out <- cudacell_seurat(
    object,
    n_hvg = 6,
    n_components = 3,
    k = 3,
    device = "cpu"
  )
  tools <- methods::slot(out, "tools")
  tools[["cudacell_seurat"]] <- list(schema = "wrong")
  methods::slot(out, "tools") <- tools
  expect_error(
    cuda_provenance(out),
    class = "cudacell_seurat_metadata_error"
  )
})

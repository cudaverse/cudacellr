.seurat_schema <- "cudacellr-seurat/1"
.seurat_tool <- "cudacell_seurat"

.seurat_stop <- function(message, class, ...) {
  fields <- list(...)
  condition <- structure(
    c(
      list(message = message, call = NULL),
      fields
    ),
    class = c(class, "cudacell_seurat_error", "error", "condition")
  )
  stop(condition)
}

.seurat_require <- function() {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    .seurat_stop(
      paste0(
        "Install the CRAN package 'SeuratObject' to use Seurat adapters. ",
        "Run `install.packages(\"SeuratObject\")`."
      ),
      class = "cudacell_seurat_dependency_error",
      missing_packages = "SeuratObject"
    )
  }
  installed <- utils::packageVersion("SeuratObject")
  required <- base::package_version("5.0.0")
  if (installed < required) {
    .seurat_stop(
      paste0(
        "Seurat adapters require SeuratObject >= 5.0.0; installed version ",
        "is ",
        as.character(installed),
        "."
      ),
      class = "cudacell_seurat_dependency_error",
      installed_version = as.character(installed),
      required_version = as.character(required)
    )
  }
  invisible(TRUE)
}

.seurat_scalar_name <- function(x, argument, allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !nzchar(x) || !identical(trimws(x), x)) {
    .seurat_stop(
      sprintf("`%s` must be one non-empty character string.", argument),
      class = "cudacell_seurat_name_error",
      argument = argument
    )
  }
  x
}

.seurat_flag <- function(x, argument) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .seurat_stop(
      sprintf("`%s` must be TRUE or FALSE.", argument),
      class = "cudacell_seurat_name_error",
      argument = argument
    )
  }
  x
}

.seurat_key <- function(name, argument) {
  value <- paste0(
    tolower(gsub("[^A-Za-z0-9]", "", name)),
    "_"
  )
  if (!grepl("^[A-Za-z][A-Za-z0-9]*_$", value)) {
    .seurat_stop(
      sprintf(
        "`%s` must contain a letter before any optional digits.",
        argument
      ),
      class = "cudacell_seurat_name_error",
      argument = argument
    )
  }
  value
}

.seurat_reduction_key <- function(key) {
  key <- .seurat_scalar_name(key, "reduction_key")
  if (!grepl("^[A-Za-z][A-Za-z0-9]*_$", key)) {
    .seurat_stop(
      paste0(
        "`reduction_key` must start with a letter, contain only letters or ",
        "digits, and end in an underscore."
      ),
      class = "cudacell_seurat_name_error",
      argument = "reduction_key"
    )
  }
  key
}

.seurat_layer_names <- function(x, assay) {
  names <- tryCatch(
    SeuratObject::Layers(x[[assay]]),
    error = function(error) {
      .seurat_stop(
        sprintf(
          "Could not inspect layers for assay '%s': %s",
          assay,
          conditionMessage(error)
        ),
        class = "cudacell_seurat_layer_error",
        assay = assay
      )
    }
  )
  as.character(names)
}

.seurat_realize_layer <- function(values, assay, layer, realize) {
  layer_class <- class(values)
  if (is.matrix(values) || methods::is(values, "Matrix")) {
    return(list(
      counts = values,
      class = layer_class,
      materialized = FALSE
    ))
  }
  if (!realize) {
    .seurat_stop(
      paste0(
        "The selected Seurat layer is not an in-memory matrix or Matrix. ",
        "cudacellr will not silently realize it; set `realize = TRUE` after ",
        "confirming that the layer fits in memory."
      ),
      class = "cudacell_seurat_layer_error",
      assay = assay,
      layer = layer,
      layer_class = layer_class
    )
  }

  original_dimnames <- tryCatch(
    dimnames(values),
    error = function(error) NULL
  )
  counts <- tryCatch(
    methods::as(values, "CsparseMatrix"),
    error = function(sparse_error) {
      tryCatch(
        Matrix::Matrix(as.matrix(values), sparse = TRUE),
        error = function(matrix_error) {
          .seurat_stop(
            sprintf(
              "Could not realize assay '%s' layer '%s' as a sparse Matrix: %s",
              assay,
              layer,
              conditionMessage(matrix_error)
            ),
            class = "cudacell_seurat_layer_error",
            assay = assay,
            layer = layer,
            layer_class = layer_class,
            sparse_error = conditionMessage(sparse_error)
          )
        }
      )
    }
  )
  if (!is.null(original_dimnames)) {
    dimnames(counts) <- original_dimnames
  }
  list(
    counts = counts,
    class = layer_class,
    materialized = TRUE
  )
}

.seurat_layer <- function(x, assay, layer, realize) {
  values <- tryCatch(
    SeuratObject::LayerData(
      x[[assay]],
      layer = layer,
      fast = FALSE
    ),
    error = function(error) {
      .seurat_stop(
        sprintf(
          "Could not read assay '%s' layer '%s': %s",
          assay,
          layer,
          conditionMessage(error)
        ),
        class = "cudacell_seurat_layer_error",
        assay = assay,
        layer = layer
      )
    }
  )
  .seurat_realize_layer(
    values,
    assay = assay,
    layer = layer,
    realize = realize
  )
}

.seurat_key_conflicts <- function(x, output_assay, reduction,
                                  assay_key, reduction_key) {
  conflicts <- character()
  subobjects <- list(
    assay = setdiff(
      SeuratObject::Assays(x),
      output_assay
    ),
    reduction = setdiff(
      SeuratObject::Reductions(x),
      reduction
    )
  )
  for (kind in names(subobjects)) {
    for (name in subobjects[[kind]]) {
      existing <- SeuratObject::Key(x[[name]])
      if (identical(existing, assay_key)) {
        conflicts <- c(
          conflicts,
          sprintf(
            "assay key '%s' used by %s '%s'",
            assay_key,
            kind,
            name
          )
        )
      }
      if (identical(existing, reduction_key)) {
        conflicts <- c(
          conflicts,
          sprintf(
            "reduction key '%s' used by %s '%s'",
            reduction_key,
            kind,
            name
          )
        )
      }
    }
  }
  conflicts
}

.seurat_cross_type_conflicts <- function(x, output_assay, reduction,
                                         neighbor) {
  all_objects <- names(x)
  expected <- list(
    output_assay = SeuratObject::Assays(x),
    reduction = SeuratObject::Reductions(x),
    neighbor = SeuratObject::Neighbors(x)
  )
  targets <- c(
    output_assay = output_assay,
    reduction = reduction,
    neighbor = neighbor
  )
  conflicts <- character()
  for (kind in names(targets)) {
    target <- targets[[kind]]
    if (target %in% all_objects && !target %in% expected[[kind]]) {
      conflicts <- c(
        conflicts,
        sprintf(
          "%s '%s' is already another Seurat subobject type",
          kind,
          target
        )
      )
    }
  }
  conflicts
}

.seurat_conflicts <- function(x, output_assay, reduction, neighbor,
                              metadata_prefix, assay_key, reduction_key) {
  conflicts <- .seurat_cross_type_conflicts(
    x,
    output_assay = output_assay,
    reduction = reduction,
    neighbor = neighbor
  )
  if (output_assay %in% SeuratObject::Assays(x)) {
    conflicts <- c(
      conflicts,
      sprintf("assay '%s'", output_assay)
    )
  }
  if (reduction %in% SeuratObject::Reductions(x)) {
    conflicts <- c(
      conflicts,
      sprintf("reduction '%s'", reduction)
    )
  }
  if (neighbor %in% SeuratObject::Neighbors(x)) {
    conflicts <- c(
      conflicts,
      sprintf("neighbor '%s'", neighbor)
    )
  }
  size_factor_field <- paste0(metadata_prefix, "_size_factor")
  if (size_factor_field %in% colnames(x[[]])) {
    conflicts <- c(
      conflicts,
      sprintf("cell metadata '%s'", size_factor_field)
    )
  }
  if (.seurat_tool %in% SeuratObject::Tool(x)) {
    conflicts <- c(
      conflicts,
      sprintf("tool '%s'", .seurat_tool)
    )
  }
  if (identical(assay_key, reduction_key)) {
    conflicts <- c(
      conflicts,
      sprintf(
        "assay key '%s' equals the requested reduction key",
        assay_key
      )
    )
  }
  c(
    conflicts,
    .seurat_key_conflicts(
      x,
      output_assay = output_assay,
      reduction = reduction,
      assay_key = assay_key,
      reduction_key = reduction_key
    )
  )
}

.seurat_hvg_metadata <- function(variable_features, feature_names,
                                 prefix) {
  required <- c(
    "index",
    "mean",
    "variance",
    "dispersion",
    "selected"
  )
  n_features <- length(feature_names)
  if (!is.data.frame(variable_features) ||
      !all(required %in% names(variable_features))) {
    .seurat_stop(
      "The workflow returned incomplete feature statistics.",
      class = "cudacell_seurat_mapping_error"
    )
  }
  positions <- variable_features$index
  if (length(positions) != n_features ||
      !is.numeric(positions) ||
      anyNA(positions) ||
      any(!is.finite(positions)) ||
      any(positions != as.integer(positions)) ||
      !identical(sort(as.integer(positions)), seq_len(n_features))) {
    .seurat_stop(
      "The workflow returned an invalid feature-position mapping.",
      class = "cudacell_seurat_mapping_error"
    )
  }
  metadata <- data.frame(
    mean = numeric(n_features),
    variance = numeric(n_features),
    dispersion = numeric(n_features),
    hvg = logical(n_features),
    hvg_rank = integer(n_features),
    row.names = feature_names,
    check.names = FALSE
  )
  metadata$mean[positions] <- variable_features$mean
  metadata$variance[positions] <- variable_features$variance
  metadata$dispersion[positions] <- variable_features$dispersion
  metadata$hvg[positions] <- variable_features$selected
  metadata$hvg_rank[positions] <- seq_len(n_features)
  names(metadata) <- paste0(prefix, "_", names(metadata))
  metadata
}

.seurat_clean_matrix <- function(x) {
  for (name in c(
    "provenance_schema",
    "compute_device",
    "compute_stages"
  )) {
    attr(x, name) <- NULL
  }
  x
}

.seurat_output_assay <- function(fit, feature_names, prefix, assay_key) {
  normalized <- .seurat_clean_matrix(fit$normalized)
  assay <- SeuratObject::CreateAssay5Object(data = normalized)
  SeuratObject::Key(assay) <- assay_key
  feature_metadata <- .seurat_hvg_metadata(
    fit$variable_features,
    feature_names = feature_names,
    prefix = prefix
  )
  assay <- SeuratObject::AddMetaData(
    assay,
    metadata = feature_metadata
  )
  selected <- fit$variable_features$index[
    fit$variable_features$selected
  ]
  SeuratObject::VariableFeatures(assay) <- feature_names[selected]
  assay
}

.seurat_reduction <- function(fit, output_assay, reduction_key) {
  SeuratObject::CreateDimReducObject(
    embeddings = fit$pca$x,
    loadings = fit$pca$rotation,
    stdev = fit$pca$sdev,
    assay = output_assay,
    key = reduction_key,
    misc = list(
      source = "cudacellr",
      provenance_schema = fit$provenance_schema,
      compute_device = fit$compute_device
    )
  )
}

.seurat_neighbor <- function(fit, cell_names) {
  index <- fit$neighbors$index
  distance <- fit$neighbors$distance
  if (!is.matrix(index) || !is.numeric(index) ||
      !is.matrix(distance) || !is.numeric(distance) ||
      !identical(dim(index), dim(distance)) ||
      nrow(index) != length(cell_names) ||
      anyNA(index) || any(!is.finite(index)) ||
      any(index != as.integer(index)) ||
      any(index < 1L | index > length(cell_names)) ||
      anyNA(distance) || any(!is.finite(distance)) ||
      any(distance < 0)) {
    .seurat_stop(
      "The workflow returned an invalid nearest-neighbour mapping.",
      class = "cudacell_seurat_mapping_error"
    )
  }
  storage.mode(index) <- "integer"
  storage.mode(distance) <- "double"
  rownames(index) <- cell_names
  rownames(distance) <- cell_names
  methods::new(
    "Neighbor",
    nn.idx = index,
    nn.dist = distance,
    alg.idx = NULL,
    alg.info = list(
      algorithm = "cudacellr_exact_knn",
      metric = fit$neighbors$metric
    ),
    cell.names = cell_names
  )
}

.seurat_set_tool <- function(x, record) {
  tools <- methods::slot(x, "tools")
  tools[[.seurat_tool]] <- record
  methods::slot(x, "tools") <- tools
  x
}

.seurat_provenance <- function(x) {
  .seurat_require()
  record <- SeuratObject::Tool(x, slot = .seurat_tool)
  if (is.null(record)) {
    .seurat_stop(
      paste0(
        "This Seurat object has no cudacellr tool record. ",
        "Run `cudacell_seurat()` first."
      ),
      class = "cudacell_seurat_metadata_error"
    )
  }
  if (!is.list(record) ||
      !identical(record$schema, .seurat_schema) ||
      is.null(record$provenance_schema) ||
      is.null(record$compute_device) ||
      is.null(record$compute_stages)) {
    .seurat_stop(
      "The stored cudacellr Seurat tool record is invalid.",
      class = "cudacell_seurat_metadata_error"
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
#' @method cuda_provenance Seurat
#' @export
cuda_provenance.Seurat <- function(x) {
  .seurat_provenance(x)
}

#' Run a cudacell workflow on a Seurat object
#'
#' `cudacell_seurat()` reads one exact assay layer, runs
#' [cudacell_workflow()], and returns a modified copy of the input object.
#' The full Seurat package is not required; the object contract uses
#' `SeuratObject` version 5 or newer.
#'
#' Results are stored in native, namespaced locations:
#'
#' - normalized expression and feature statistics in a new `Assay5`;
#' - PCA embeddings and loadings in a `DimReduc`;
#' - exact kNN indices and distances in a `Neighbor`;
#' - library-size factors in cell metadata;
#' - parameters and compute provenance in the `cudacell_seurat` tool record.
#'
#' Existing assays and layers, reductions, graphs, neighbours, identities,
#' cell metadata, images, commands, miscellaneous data, and tool records are
#' retained. Every output collision is checked before device selection or
#' computation. `overwrite = TRUE` replaces only the explicitly named
#' cudacellr outputs.
#'
#' Feature metadata use the requested prefix. `<prefix>_hvg_rank` is the
#' one-based global dispersion rank across all input features;
#' `<prefix>_hvg` identifies the selected top `n_hvg` features.
#'
#' Non-memory-backed layers are never realized silently. Set `realize = TRUE`
#' only after confirming that the selected layer fits in memory. The selected
#' layer must contain every object cell in the same order.
#'
#' @param x A `Seurat` object.
#' @param assay Exact input assay name. `NULL` uses
#'   [SeuratObject::DefaultAssay()].
#' @param layer Exact feature-by-cell count layer name.
#' @param n_hvg Number of highly variable features.
#' @param n_components Number of PCA components.
#' @param k Number of neighbours per cell.
#' @param device Computation device.
#' @param batch_size Maximum query rows per exact kNN distance block.
#' @param scale_factor Target library size used during normalization.
#' @param log1p Whether to apply `log1p()` after library-size normalization.
#' @param min_mean Minimum feature mean used for HVG selection.
#' @param scale. Whether to scale selected features before PCA.
#' @param output_assay Name for the added normalized `Assay5`. Its `data`
#'   layer contains normalized expression.
#' @param reduction Name for the added PCA `DimReduc`.
#' @param reduction_key Seurat key for PCA dimensions. It must start with a
#'   letter, contain only letters or digits, and end in an underscore.
#' @param neighbor Name for the added exact kNN `Neighbor`.
#' @param metadata_prefix Prefix for feature statistics and the cell-level
#'   size-factor column.
#' @param overwrite Whether the named cudacellr output fields may be replaced.
#'   Other object state is never overwritten.
#' @param realize Whether a non-memory-backed layer may be explicitly
#'   materialized as a sparse `Matrix`.
#' @return A valid `Seurat` object. The input object is not modified.
#' @export
#' @examples
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   counts <- Matrix::Matrix(
#'     matrix(
#'       rpois(30 * 12, lambda = 2),
#'       nrow = 30,
#'       dimnames = list(
#'         paste0("gene_", 1:30),
#'         paste0("cell_", 1:12)
#'       )
#'     ),
#'     sparse = TRUE
#'   )
#'   object <- SeuratObject::CreateSeuratObject(counts)
#'   object <- cudacell_seurat(
#'     object,
#'     n_hvg = 10,
#'     n_components = 3,
#'     k = 3,
#'     device = "cpu"
#'   )
#'   SeuratObject::Embeddings(object[["cudacell_pca"]])
#'   cuda_provenance(object)
#' }
cudacell_seurat <- function(x, assay = NULL, layer = "counts",
                            n_hvg = 2000L, n_components = 30L, k = 15L,
                            device = c("auto", "cuda", "cpu"),
                            batch_size = 256L, scale_factor = 10000,
                            log1p = TRUE, min_mean = 0, scale. = TRUE,
                            output_assay = "CUDACELL",
                            reduction = "cudacell_pca",
                            reduction_key = "CUDACELLPC_",
                            neighbor = "cudacell_knn",
                            metadata_prefix = "cudacell",
                            overwrite = FALSE, realize = FALSE) {
  .seurat_require()
  if (!inherits(x, "Seurat")) {
    .seurat_stop(
      "`x` must be a Seurat object.",
      class = "cudacell_seurat_type_error"
    )
  }
  assay <- .seurat_scalar_name(assay, "assay", allow_null = TRUE)
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(x)
  }
  assay <- .seurat_scalar_name(assay, "assay")
  layer <- .seurat_scalar_name(layer, "layer")
  output_assay <- .seurat_scalar_name(
    output_assay,
    "output_assay"
  )
  reduction <- .seurat_scalar_name(reduction, "reduction")
  neighbor <- .seurat_scalar_name(neighbor, "neighbor")
  metadata_prefix <- .seurat_scalar_name(
    metadata_prefix,
    "metadata_prefix"
  )
  reduction_key <- .seurat_reduction_key(reduction_key)
  log1p <- .seurat_flag(log1p, "log1p")
  scale. <- .seurat_flag(scale., "scale.")
  overwrite <- .seurat_flag(overwrite, "overwrite")
  realize <- .seurat_flag(realize, "realize")

  assay_names <- SeuratObject::Assays(x)
  if (!assay %in% assay_names) {
    .seurat_stop(
      sprintf(
        "`assay` '%s' was not found. Available assays: %s.",
        assay,
        paste(sprintf("'%s'", assay_names), collapse = ", ")
      ),
      class = "cudacell_seurat_layer_error",
      assay = assay,
      available_assays = assay_names
    )
  }
  layer_names <- .seurat_layer_names(x, assay)
  if (!layer %in% layer_names) {
    .seurat_stop(
      sprintf(
        "`layer` '%s' was not found in assay '%s'. Available layers: %s.",
        layer,
        assay,
        paste(sprintf("'%s'", layer_names), collapse = ", ")
      ),
      class = "cudacell_seurat_layer_error",
      assay = assay,
      layer = layer,
      available_layers = layer_names
    )
  }
  if (identical(output_assay, assay)) {
    .seurat_stop(
      "`output_assay` must differ from the input `assay`.",
      class = "cudacell_seurat_name_error",
      argument = "output_assay"
    )
  }
  output_names <- c(output_assay, reduction, neighbor)
  if (anyDuplicated(output_names)) {
    .seurat_stop(
      "`output_assay`, `reduction`, and `neighbor` must use distinct names.",
      class = "cudacell_seurat_name_error"
    )
  }
  assay_key <- .seurat_key(output_assay, "output_assay")
  conflicts <- .seurat_conflicts(
    x,
    output_assay = output_assay,
    reduction = reduction,
    neighbor = neighbor,
    metadata_prefix = metadata_prefix,
    assay_key = assay_key,
    reduction_key = reduction_key
  )
  permanent_conflicts <- grep(
    "another Seurat subobject type| key '",
    conflicts,
    value = TRUE
  )
  if (length(permanent_conflicts)) {
    .seurat_stop(
      paste0(
        "The requested output names or keys conflict with unrelated Seurat ",
        "state: ",
        paste(permanent_conflicts, collapse = ", "),
        "."
      ),
      class = "cudacell_seurat_collision_error",
      conflicts = permanent_conflicts
    )
  }
  if (length(conflicts) && !overwrite) {
    .seurat_stop(
      paste0(
        "cudacellr output fields already exist: ",
        paste(conflicts, collapse = ", "),
        ". Set `overwrite = TRUE` to replace only these fields."
      ),
      class = "cudacell_seurat_collision_error",
      conflicts = conflicts
    )
  }

  device <- match.arg(device)
  cudaverse::cuda_select_device(device)
  input <- .seurat_layer(
    x,
    assay = assay,
    layer = layer,
    realize = realize
  )
  counts <- .cell_counts(input$counts)
  feature_names <- rownames(counts)
  cell_names <- colnames(counts)
  if (is.null(feature_names) || anyDuplicated(feature_names)) {
    .seurat_stop(
      "The selected layer must have unique feature names.",
      class = "cudacell_seurat_mapping_error"
    )
  }
  if (is.null(cell_names) || anyDuplicated(cell_names) ||
      !identical(cell_names, colnames(x))) {
    .seurat_stop(
      paste0(
        "The selected layer must contain every Seurat object cell exactly ",
        "once and in object order."
      ),
      class = "cudacell_seurat_mapping_error"
    )
  }

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
  if (!identical(dim(fit$normalized), dim(counts)) ||
      !identical(dimnames(fit$normalized), dimnames(counts)) ||
      nrow(fit$pca$x) != length(cell_names) ||
      !identical(rownames(fit$pca$x), cell_names) ||
      nrow(fit$pca$rotation) != sum(fit$variable_features$selected)) {
    .seurat_stop(
      "The workflow results do not align with the Seurat input layer.",
      class = "cudacell_seurat_mapping_error"
    )
  }

  size_factors <- as.numeric(Matrix::colSums(counts)) / scale_factor
  names(size_factors) <- cell_names
  size_factor_field <- paste0(metadata_prefix, "_size_factor")
  cell_metadata <- data.frame(
    size_factors,
    row.names = cell_names,
    check.names = FALSE
  )
  names(cell_metadata) <- size_factor_field

  out <- x
  output_assay_object <- .seurat_output_assay(
    fit,
    feature_names = feature_names,
    prefix = metadata_prefix,
    assay_key = assay_key
  )
  out <- withCallingHandlers(
    {
      out[[output_assay]] <- output_assay_object
      out
    },
    warning = function(condition) {
      if (identical(
        conditionMessage(condition),
        "No layers found matching search pattern provided"
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
  out <- SeuratObject::AddMetaData(
    out,
    metadata = cell_metadata
  )
  reduction_object <- .seurat_reduction(
    fit,
    output_assay = output_assay,
    reduction_key = reduction_key
  )
  out <- withCallingHandlers(
    {
      out[[reduction]] <- reduction_object
      out
    },
    warning = function(condition) {
      if (grepl(
        "^Number of dimensions changing from [0-9]+ to [0-9]+$",
        conditionMessage(condition)
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
  out[[neighbor]] <- .seurat_neighbor(
    fit,
    cell_names = cell_names
  )
  record <- list(
    schema = .seurat_schema,
    package_version = as.character(
      utils::packageVersion("cudacellr")
    ),
    input_assay = assay,
    input_layer = layer,
    input_layer_class = input$class,
    input_layer_materialized = input$materialized,
    outputs = list(
      assay = output_assay,
      layer = "data",
      reduction = reduction,
      reduction_key = reduction_key,
      neighbor = neighbor,
      feature_metadata = paste0(
        metadata_prefix,
        "_",
        c("mean", "variance", "dispersion", "hvg", "hvg_rank")
      ),
      size_factor = size_factor_field
    ),
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
  out <- .seurat_set_tool(out, record)

  methods::validObject(out)
  out
}

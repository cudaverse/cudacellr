# Run a cudacell workflow on a Seurat object

`cudacell_seurat()` reads one exact assay layer, runs
[`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md),
and returns a modified copy of the input object. The full Seurat package
is not required; the object contract uses `SeuratObject` version 5 or
newer.

## Usage

``` r
cudacell_seurat(
  x,
  assay = NULL,
  layer = "counts",
  n_hvg = 2000L,
  n_components = 30L,
  k = 15L,
  device = c("auto", "cuda", "cpu"),
  batch_size = 256L,
  scale_factor = 10000,
  log1p = TRUE,
  min_mean = 0,
  scale. = TRUE,
  output_assay = "CUDACELL",
  reduction = "cudacell_pca",
  reduction_key = "CUDACELLPC_",
  neighbor = "cudacell_knn",
  metadata_prefix = "cudacell",
  overwrite = FALSE,
  realize = FALSE
)
```

## Arguments

- x:

  A `Seurat` object.

- assay:

  Exact input assay name. `NULL` uses
  [`SeuratObject::DefaultAssay()`](https://satijalab.github.io/seurat-object/reference/DefaultAssay.html).

- layer:

  Exact feature-by-cell count layer name.

- n_hvg:

  Number of highly variable features.

- n_components:

  Number of PCA components.

- k:

  Number of neighbours per cell.

- device:

  Computation device.

- batch_size:

  Maximum query rows per exact kNN distance block.

- scale_factor:

  Target library size used during normalization.

- log1p:

  Whether to apply [`log1p()`](https://rdrr.io/r/base/Log.html) after
  library-size normalization.

- min_mean:

  Minimum feature mean used for HVG selection.

- scale.:

  Whether to scale selected features before PCA.

- output_assay:

  Name for the added normalized `Assay5`. Its `data` layer contains
  normalized expression.

- reduction:

  Name for the added PCA `DimReduc`.

- reduction_key:

  Seurat key for PCA dimensions. It must start with a letter, contain
  only letters or digits, and end in an underscore.

- neighbor:

  Name for the added exact kNN `Neighbor`.

- metadata_prefix:

  Prefix for feature statistics and the cell-level size-factor column.

- overwrite:

  Whether the named cudacellr output fields may be replaced. Other
  object state is never overwritten.

- realize:

  Whether a non-memory-backed layer may be explicitly materialized as a
  sparse `Matrix`.

## Value

A valid `Seurat` object. The input object is not modified.

## Details

Results are stored in native, namespaced locations:

- normalized expression and feature statistics in a new `Assay5`;

- PCA embeddings and loadings in a `DimReduc`;

- exact kNN indices and distances in a `Neighbor`;

- library-size factors in cell metadata;

- parameters and compute provenance in the `cudacell_seurat` tool
  record.

Existing assays and layers, reductions, graphs, neighbours, identities,
cell metadata, images, commands, miscellaneous data, and tool records
are retained. Every output collision is checked before device selection
or computation. `overwrite = TRUE` replaces only the explicitly named
cudacellr outputs.

Feature metadata use the requested prefix. `<prefix>_hvg_rank` is the
one-based global dispersion rank across all input features;
`<prefix>_hvg` identifies the selected top `n_hvg` features.

Non-memory-backed layers are never realized silently. Set
`realize = TRUE` only after confirming that the selected layer fits in
memory. The selected layer must contain every object cell in the same
order.

## Examples

``` r
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  counts <- Matrix::Matrix(
    matrix(
      rpois(30 * 12, lambda = 2),
      nrow = 30,
      dimnames = list(
        paste0("gene_", 1:30),
        paste0("cell_", 1:12)
      )
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts)
  object <- cudacell_seurat(
    object,
    n_hvg = 10,
    n_components = 3,
    k = 3,
    device = "cpu"
  )
  SeuratObject::Embeddings(object[["cudacell_pca"]])
  cuda_provenance(object)
}
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> <cuda_provenance schema=cudaverse-stage/1 stages=6 compute=cpu>
#>                   stage requested_device device backend   selection_reason
#>           normalization        fixed-cpu    cpu  Matrix algorithm_cpu_only
#>                     hvg        fixed-cpu    cpu  Matrix algorithm_cpu_only
#>       pca_preprocessing              cpu    cpu   stats       explicit_cpu
#>       pca_decomposition              cpu    cpu   stats       explicit_cpu
#>            knn_distance              cpu    cpu    base       explicit_cpu
#>  knn_neighbor_selection        fixed-cpu    cpu    base algorithm_cpu_only
#>  fallback output_device
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
```

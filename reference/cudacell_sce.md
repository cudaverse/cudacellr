# Run a cudacell workflow on a SingleCellExperiment

`cudacell_sce()` reads one assay, runs
[`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md),
and returns a modified copy of the input object. Existing assays, row
and column metadata, reduced dimensions, alternative experiments,
pairings, and top-level metadata are retained. New results are stored in
native `SingleCellExperiment` locations:

## Usage

``` r
cudacell_sce(
  x,
  assay = "counts",
  n_hvg = 2000L,
  n_components = 30L,
  k = 15L,
  device = c("auto", "cuda", "cpu"),
  batch_size = 256L,
  scale_factor = 10000,
  log1p = TRUE,
  min_mean = 0,
  scale. = TRUE,
  normalized_assay = NULL,
  reduced_dim = "CUDACELL_PCA",
  col_pair = "CUDACELL_KNN",
  row_data_prefix = "cudacell",
  set_size_factors = FALSE,
  overwrite = FALSE,
  realize = FALSE
)
```

## Arguments

- x:

  A `SingleCellExperiment`.

- assay:

  Exact name of the feature-by-cell count assay.

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

- normalized_assay:

  Name for the added normalized assay. `NULL` chooses
  `"cudacell_logcounts"` when `log1p = TRUE` and `"cudacell_normalized"`
  otherwise.

- reduced_dim:

  Name for the added PCA reduced dimension.

- col_pair:

  Name for the added directed kNN `colPair`.

- row_data_prefix:

  Prefix for added feature statistics and the namespaced column size
  factor.

- set_size_factors:

  Whether to also replace the canonical
  [`SingleCellExperiment::sizeFactors()`](https://rdrr.io/pkg/BiocGenerics/man/dge.html)
  values. The namespaced `<row_data_prefix>_size_factor` column is
  always added to `colData`.

- overwrite:

  Whether explicitly named cudacellr output fields may be replaced.
  Other object content is never overwritten.

- realize:

  Whether a delayed assay may be explicitly materialized as a sparse
  `Matrix`.

## Value

A valid `SingleCellExperiment` of the same class as `x`.

## Details

- normalized expression in an assay;

- PCA scores in a reduced dimension;

- highly variable feature statistics in `rowData`;

- directed nearest-neighbour relationships in a `colPair`;

- parameters and compute provenance in `metadata(x)$cudacellr`.

The selected assay must currently be a base matrix or `Matrix` object.
Delayed assays are never realized silently. Set `realize = TRUE` only
after confirming that the selected assay fits in memory.

## Examples

``` r
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  counts <- matrix(
    rpois(30 * 12, lambda = 2),
    nrow = 30,
    dimnames = list(paste0("gene_", 1:30), paste0("cell_", 1:12))
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )
  sce <- cudacell_sce(
    sce,
    n_hvg = 10,
    n_components = 3,
    k = 3,
    device = "cpu"
  )
  SingleCellExperiment::reducedDim(sce, "CUDACELL_PCA")
  cuda_provenance(sce)
}
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

# A native SeuratObject v5 workflow

[`cudacell_seurat()`](https://cudaverse.github.io/cudacellr/reference/cudacell_seurat.md)
runs the complete cudacellr preprocessing, PCA, and exact
nearest-neighbour workflow while keeping every result in a native Seurat
v5 subobject. The adapter requires only `SeuratObject >= 5.0.0`;
installing the full Seurat package is not necessary.

``` r

install.packages("SeuratObject")
```

## Start with a populated object

This example contains a second assay, cell metadata, identities, a
reduction, miscellaneous data, and a tool record. The workflow returns a
modified copy and retains all of them.

``` r

library(cudacellr)
library(Matrix)

set.seed(42)
counts <- Matrix(
  matrix(
    rpois(80 * 30, lambda = 2),
    nrow = 80,
    dimnames = list(
      paste0("gene_", seq_len(80)),
      paste0("cell_", seq_len(30))
    )
  ),
  sparse = TRUE
)
counts[1, ] <- counts[1, ] + 1

object <- SeuratObject::CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  meta.data = data.frame(
    batch = rep(c("one", "two"), each = 15),
    row.names = colnames(counts)
  ),
  project = "preserved-project"
)
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes
#> ('-')
object[["ADT"]] <- SeuratObject::CreateAssay5Object(
  counts = counts[seq_len(5), , drop = FALSE]
)
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes
#> ('-')
object[["existing_pca"]] <- SeuratObject::CreateDimReducObject(
  embeddings = matrix(
    rnorm(ncol(counts) * 2),
    nrow = ncol(counts),
    dimnames = list(colnames(counts), c("OLDPC_1", "OLDPC_2"))
  ),
  assay = "RNA",
  key = "OLDPC_"
)
SeuratObject::Idents(object) <- factor(
  rep(c("A", "B"), each = 15)
)
SeuratObject::Misc(object, slot = "owner") <- "unchanged"
add_example_tool <- function(x) {
  SeuratObject::Tool(x) <- list(owner = "unchanged")
  x
}
object <- add_example_tool(object)
existing_tool_names <- SeuratObject::Tool(object)
```

## Run the workflow

Input assay and layer names are exact. `assay = NULL` uses the default
assay, while `layer = "counts"` is explicit by default.

``` r

result <- cudacell_seurat(
  object,
  assay = "RNA",
  layer = "counts",
  n_hvg = 30,
  n_components = 5,
  k = 5,
  batch_size = 8,
  device = "cpu"
)

SeuratObject::Assays(result)
#> [1] "RNA"      "ADT"      "CUDACELL"
SeuratObject::Reductions(result)
#> [1] "existing_pca" "cudacell_pca"
SeuratObject::Neighbors(result)
#> [1] "cudacell_knn"
```

The default outputs use namespaced native locations:

``` r

normalized <- SeuratObject::LayerData(
  result[["CUDACELL"]],
  layer = "data",
  fast = FALSE
)
normalized[1:3, 1:3]
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>          cell_1   cell_2   cell_3
#> gene-1 5.652808 5.208492 5.773037
#> gene-2 5.430541 4.120760 4.861401
#> gene-3 4.057303 4.120760 4.861401

SeuratObject::Embeddings(result[["cudacell_pca"]])[1:3, , drop = FALSE]
#>        CUDACELLPC_1 CUDACELLPC_2 CUDACELLPC_3 CUDACELLPC_4 CUDACELLPC_5
#> cell_1    -1.010505     2.472902    -1.934263    1.0690482  -0.05279964
#> cell_2    -1.130595     1.670652     1.065843   -2.2639677  -3.15504850
#> cell_3     2.060318     3.225506    -1.481663    0.6235578  -0.62236913
SeuratObject::Indices(result[["cudacell_knn"]])[1:3, , drop = FALSE]
#>        neighbor_1 neighbor_2 neighbor_3 neighbor_4 neighbor_5
#> cell_1          9         18         17          8         29
#> cell_2          6         21         29         28         23
#> cell_3         18         29         26          1          9
result[["CUDACELL"]][[]][
  1:3,
  c(
    "cudacell_mean",
    "cudacell_variance",
    "cudacell_dispersion",
    "cudacell_hvg",
    "cudacell_hvg_rank"
  )
]
#>        cudacell_mean cudacell_variance cudacell_dispersion cudacell_hvg
#> gene-1      5.241438         0.2898075           0.0552916        FALSE
#> gene-2      4.517405         1.7976247           0.3979330        FALSE
#> gene-3      4.028052         2.8141622           0.6986410        FALSE
#>        cudacell_hvg_rank
#> gene-1                80
#> gene-2                73
#> gene-3                50
```

`cudacell_hvg_rank` is the one-based global dispersion rank across every
input feature. `cudacell_hvg` marks the selected top `n_hvg` features,
which are also registered through
[`SeuratObject::VariableFeatures()`](https://satijalab.github.io/seurat-object/reference/VariableFeatures.html).

The exact directed kNN result is stored as a `Neighbor`, including both
integer indices and distances. The PCA `DimReduc` includes cell
embeddings, feature loadings, standard deviations, and the `CUDACELLPC_`
key.

## Verify preservation

The input object is unchanged, and unrelated state survives in the
result:

``` r

identical(
  SeuratObject::LayerData(object[["RNA"]], "counts"),
  SeuratObject::LayerData(result[["RNA"]], "counts")
)
#> [1] TRUE
identical(object[["ADT"]], result[["ADT"]])
#> [1] TRUE
identical(object[["existing_pca"]], result[["existing_pca"]])
#> [1] TRUE
identical(SeuratObject::Idents(object), SeuratObject::Idents(result))
#> [1] TRUE
identical(
  SeuratObject::Misc(object, slot = "owner"),
  SeuratObject::Misc(result, slot = "owner")
)
#> [1] TRUE
all(vapply(
  existing_tool_names,
  function(name) {
    identical(
      SeuratObject::Tool(object, slot = name),
      SeuratObject::Tool(result, slot = name)
    )
  },
  logical(1)
))
#> [1] TRUE
```

cudacellr adds its own deterministic `cudacell_seurat` tool record
without replacing existing tool records.

## Collision and realization safety

Every target assay, reduction, neighbour, feature-metadata column,
cell-metadata column, tool record, and Seurat key is checked before
device selection, layer access, realization, or numerical work. A second
run with the same names fails:

``` r

cudacell_seurat(
  result,
  n_hvg = 30,
  n_components = 5,
  k = 5,
  device = "cpu"
)
#> Error:
#> ! cudacellr output fields already exist: assay 'CUDACELL', reduction 'cudacell_pca', neighbor 'cudacell_knn', cell metadata 'cudacell_size_factor', tool 'cudacell_seurat'. Set `overwrite = TRUE` to replace only these fields.
```

Use `overwrite = TRUE` only when replacing those exact cudacellr outputs
is intentional. Cross-type name or key collisions are never overwritten.

Only ordinary in-memory `matrix` and `Matrix` layers are accepted by
default. Delayed, on-disk, and other matrix-like layers are not silently
collected. Set `realize = TRUE` only after confirming that the selected
layer fits in memory; the original layer class and materialization
decision are recorded.

For `device = "cuda"`, CUDA availability is validated before the layer
is read or realized. The request remains strict and never silently falls
back to CPU.

## Provenance

The native tool record contains identifiers, parameters, output
locations, and the same ordered `cudaverse-stage/1` compute provenance
used by the matrix and SingleCellExperiment workflows:

``` r

cuda_provenance(result)
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
SeuratObject::Tool(result, slot = "cudacell_seurat")$outputs
#> $assay
#> [1] "CUDACELL"
#> 
#> $layer
#> [1] "data"
#> 
#> $reduction
#> [1] "cudacell_pca"
#> 
#> $reduction_key
#> [1] "CUDACELLPC_"
#> 
#> $neighbor
#> [1] "cudacell_knn"
#> 
#> $feature_metadata
#> [1] "cudacell_mean"       "cudacell_variance"   "cudacell_dispersion"
#> [4] "cudacell_hvg"        "cudacell_hvg_rank"  
#> 
#> $size_factor
#> [1] "cudacell_size_factor"
SeuratObject::Tool(result, slot = "cudacell_seurat")$parameters
#> $n_hvg
#> [1] 30
#> 
#> $n_components
#> [1] 5
#> 
#> $k
#> [1] 5
#> 
#> $requested_device
#> [1] "cpu"
#> 
#> $batch_size
#> [1] 8
#> 
#> $scale_factor
#> [1] 10000
#> 
#> $log1p
#> [1] TRUE
#> 
#> $min_mean
#> [1] 0
#> 
#> $scale
#> [1] TRUE
```

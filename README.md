# cudacellr

`cudacellr` is the R-native single-cell workflow layer of the **cudaverse**.
It composes sparse count processing with GPU-aware numerical algorithms while
keeping features in rows and cells in columns at the public R boundary.

## Current workflow

- library-size normalization;
- sparse log-normalization;
- highly variable feature selection;
- PCA through `cudalearnr`;
- k-nearest neighbours through `cudalearnr`;
- a composable end-to-end workflow result that reuses each preprocessing
  stage instead of repeating normalization and feature selection.

## Installation

```r
# install.packages("pak")
pak::pak("cudaverse/cudacellr")
```

## Example

```r
library(cudacellr)
library(Matrix)

set.seed(1)
counts <- Matrix(
  matrix(rpois(1000 * 300, lambda = 1.5), 1000, 300),
  sparse = TRUE
)

fit <- cudacell_workflow(
  counts,
  n_hvg = 300,
  n_components = 20,
  k = 15,
  batch_size = 128
)

fit$pca
fit$neighbors
```

Feature and cell identifiers are preserved throughout the workflow:

```r
identical(rownames(fit$normalized), rownames(counts))
identical(colnames(fit$normalized), colnames(counts))
identical(rownames(fit$pca$x), colnames(counts))
identical(rownames(fit$neighbors$index), colnames(counts))
```

HVG results also include the original feature `index`, so repeated feature
names do not cause PCA to select the wrong matrix row.

Continue from this result through graph clustering and embedding in the
cudaverse
[end-to-end workflow](https://github.com/cudaverse/.github/blob/main/WORKFLOW.md).

Neighbour search is exact but uses bounded distance blocks. At most
`min(batch_size, cells) * cells` distances are held at once instead of a full
cell-by-cell matrix. Reduce `batch_size` when memory is constrained; the
selected neighbours are deterministic and do not change with batch size.

## Backend provenance

Normalization and highly variable feature selection stay sparse and run on
the CPU. PCA and kNN distance blocks can use CUDA, while deterministic
neighbour selection remains on the CPU:

| Function | Device-selected work | Always-CPU work | CUDA aggregate |
|---|---|---|---|
| `cuda_normalize_counts()` | none | sparse normalization | `cpu` |
| `cuda_hvg()` | none | sparse feature statistics and ranking | `cpu` |
| `cuda_cell_pca()` | PCA preprocessing and decomposition | normalization and HVG selection | `hybrid` |
| `cuda_cell_neighbors()` | kNN distance blocks | neighbour selection | `hybrid` |
| `cudacell_workflow()` | PCA and kNN distance stages | sparse preprocessing and neighbour selection | `hybrid` |

Inspect the runtime and any result without guessing from its function name:

```r
cudatensr::cuda_diagnostics()
cuda_provenance(fit)
```

The provenance table separates the requested device, actual compute device,
output device, and an automatic fallback. See
[Backend provenance and CUDA diagnostics](https://cudaverse.github.io/cudacellr/articles/backend-provenance.html)
for a runnable CPU workflow, dense-PCA and kNN memory guidance, and the
hardware-CI gate.

## Object integration

The numerical core accepts base, `Matrix`, and `cudasparse` matrices.
`cudacell_sce()` adds an optional native `SingleCellExperiment` workflow
without making Bioconductor packages hard dependencies:

```r
library(SingleCellExperiment)

sce <- SingleCellExperiment(assays = list(counts = counts))
sce <- cudacell_sce(
  sce,
  n_hvg = 300,
  n_components = 20,
  k = 15,
  batch_size = 128,
  device = "cpu"
)

assayNames(sce)
reducedDimNames(sce)
colPairNames(sce)
cuda_provenance(sce)
```

The adapter writes natural-log normalized expression to
`cudacell_logcounts`, PCA scores to `CUDACELL_PCA`, HVG statistics to
namespaced `rowData` columns, and directed kNN relationships to
`CUDACELL_KNN`. Existing assays, metadata, reduced dimensions, alternative
experiments, labels, pairings, and size factors remain unchanged.

Output names are deliberately namespaced and never silently replaced.
`overwrite = TRUE` replaces only the explicitly named cudacellr fields.
Delayed assays also remain lazy unless `realize = TRUE` is requested after
checking that the assay fits in memory. See the
[SingleCellExperiment workflow](https://cudaverse.github.io/cudacellr/articles/single-cell-experiment.html)
for the complete object contract.

Seurat v5 objects have the same non-destructive workflow through
`cudacell_seurat()`. Only the lightweight `SeuratObject` package is needed;
the full Seurat package is not required:

```r
object <- SeuratObject::CreateSeuratObject(counts = counts)
object <- cudacell_seurat(
  object,
  assay = "RNA",
  layer = "counts",
  n_hvg = 300,
  n_components = 20,
  k = 15,
  batch_size = 128,
  device = "cpu"
)

SeuratObject::Assays(object)
SeuratObject::Reductions(object)
SeuratObject::Neighbors(object)
cuda_provenance(object)
```

Normalized expression and feature statistics are written to a native
`Assay5`, PCA scores and loadings to a `DimReduc`, exact neighbour indices
and distances to a `Neighbor`, size factors to cell metadata, and parameters
and compute provenance to the `cudacell_seurat` tool record. Existing assays,
layers, reductions, graphs, neighbours, identities, metadata, and tools are
preserved.

All output names and Seurat keys are checked before device selection,
realization, or computation. `overwrite = TRUE` targets only the named
cudacellr outputs, and non-memory-backed layers require the explicit
`realize = TRUE` opt-in. See the
[SeuratObject v5 workflow](https://cudaverse.github.io/cudacellr/articles/seurat-object.html)
for the full contract.

For installation, device verification, memory advice, and common failures, see
the cudaverse
[GPU setup and troubleshooting guide](https://github.com/cudaverse/.github/blob/main/GPU_SETUP.md).

## License

MIT © Yaoxiang Li

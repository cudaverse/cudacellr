# cudacellr

`cudacellr` takes a single-cell count matrix through normalization, variable
feature selection, PCA, and exact nearest neighbours. Use it to accelerate
PCA and neighbour search with NVIDIA CUDA while keeping familiar R matrices,
cell identifiers, and optional SingleCellExperiment or SeuratObject results.

The native CUDA backend comes from **cudaverse** and does not require torch or
LibTorch. Sparse normalization and feature selection run on the CPU; PCA,
distance calculations, and deterministic top-k selection can run on CUDA.

## Current workflow

- library-size normalization;
- sparse log-normalization;
- highly variable feature selection;
- PCA through `cudaverse`;
- k-nearest neighbours through `cudaverse`;
- a composable end-to-end workflow result that reuses each preprocessing
  stage instead of repeating normalization and feature selection.

## Installation

For CUDA execution, prepare a Windows or Linux machine with an NVIDIA driver,
cuBLAS 12, and cuSOLVER 11. Follow the
[CUDA setup guide](https://cudaverse.github.io/cudaverse/articles/gpu-setup.html)
before running the GPU examples. Current CUDA execution is unavailable on
macOS; use a Windows or Linux GPU machine for these workflows.

```r
# install.packages("pak")
pak::pak("cudaverse/cudacellr")
library(cudacellr)

health <- cudaverse::cuda_diagnostics()
health$summary
health$next_steps
health$selected_backend  # "native" for the lightweight CUDA path
cudaverse::cuda_select_device("cuda")
```

An explicit `device = "cuda"` request fails with diagnostics when CUDA is
unavailable; it does not silently run the requested GPU stages on the CPU.

## A single-cell CUDA workflow

Provide features in rows and cells in columns:

```r
library(cudacellr)
library(Matrix)

set.seed(1)
counts <- Matrix(
  matrix(
    rpois(1000 * 300, lambda = 1.5), 1000, 300,
    dimnames = list(
      paste0("gene", seq_len(1000)),
      paste0("cell", seq_len(300))
    )
  ),
  sparse = TRUE
)

fit <- cudacell_workflow(
  counts,
  n_hvg = 300,
  n_components = 20,
  k = 15,
  batch_size = 128,
  device = "cuda"
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

With the native backend, the compute boundaries are:

| Function | Native CUDA work | CPU work and R outputs |
|---|---|---|
| `cuda_normalize_counts()` | none | sparse normalization |
| `cuda_hvg()` | none | sparse feature statistics and ranking |
| `cuda_cell_pca()` | PCA centring/scaling and decomposition | normalization, HVG selection, selected dense input, and the public PCA model |
| `cuda_cell_neighbors()` | distance blocks and deterministic top-k selection | input validation and final neighbour matrices |
| `cudacell_workflow()` | PCA and exact kNN | sparse preprocessing and public R results |

The complete workflow is `hybrid` because normalization and HVG selection are
CPU stages. Native PCA returns ordinary R scores and loadings and also retains
a device-side score cache; passing the PCA result directly to
`cuda_cell_neighbors()` lets kNN reuse that cache without uploading the scores
again. CUDA neighbour selection is not a fixed CPU stage.

The optional torch compatibility backend can have different stage boundaries.
Check provenance for the backend actually selected rather than inferring it
from the requested device.

Inspect the runtime and any result without guessing from its function name:

```r
cudaverse::cuda_diagnostics()
cuda_provenance(fit)
```

The provenance table separates the requested device, actual compute device,
output device, and an automatic fallback. See
[Backend provenance and CUDA diagnostics](https://cudaverse.github.io/cudacellr/articles/backend-provenance.html)
for a worked CUDA example, a small CPU reference, and dense-PCA/kNN memory
guidance.

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
  device = "cuda"
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
  device = "cuda"
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

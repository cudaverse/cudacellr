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

Continue from this result through graph clustering and embedding in the
cudaverse
[end-to-end workflow](https://github.com/cudaverse/.github/blob/main/WORKFLOW.md).

Neighbour search is exact but uses bounded distance blocks. At most
`min(batch_size, cells) * cells` distances are held at once instead of a full
cell-by-cell matrix. Reduce `batch_size` when memory is constrained; the
selected neighbours are deterministic and do not change with batch size.

## Object integration

The current API accepts base, `Matrix`, and `cudasparse` matrices. Native
SingleCellExperiment and Seurat adapters are the next integration milestone;
they are intentionally not hard dependencies of the numerical core.

For installation, device verification, memory advice, and common failures, see
the cudaverse
[GPU setup and troubleshooting guide](https://github.com/cudaverse/.github/blob/main/GPU_SETUP.md).

## License

MIT © Yaoxiang Li

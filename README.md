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
  k = 15
)

fit$pca
fit$neighbors
```

## Object integration

The 0.1.0 API accepts base, `Matrix`, and `cudasparse` matrices. Native
SingleCellExperiment and Seurat adapters are the next integration milestone;
they are intentionally not hard dependencies of the numerical core.

## License

MIT © Yaoxiang Li

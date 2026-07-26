# Run the initial cudacellr workflow

Run the initial cudacellr workflow

## Usage

``` r
cudacell_workflow(
  counts,
  n_hvg = 2000L,
  n_components = 30L,
  k = 15L,
  device = c("auto", "cuda", "cpu"),
  batch_size = 256L,
  scale_factor = 10000,
  log1p = TRUE,
  min_mean = 0,
  scale. = TRUE
)
```

## Arguments

- counts:

  Feature-by-cell count matrix.

- n_hvg:

  Number of variable features.

- n_components:

  PCA components.

- k:

  Neighbours per cell.

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

  Minimum feature mean used during highly variable feature selection.

- scale.:

  Whether to scale selected features before PCA.

## Value

A list containing normalized counts, variable features, PCA, and kNN.
Feature and cell identifiers are retained throughout all stages.

## Examples

``` r
set.seed(3)
counts <- matrix(rpois(40 * 25, 2), 40, 25)
cudacell_workflow(
  counts, n_hvg = 15, n_components = 5, k = 5, device = "cpu"
)
#> <cudacell_workflow features=40 cells=25 hvg=15 components=5 k=5 pca_device=cpu compute=cpu>
```

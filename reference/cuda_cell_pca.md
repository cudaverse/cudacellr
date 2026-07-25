# Run PCA on highly variable single-cell features

Run PCA on highly variable single-cell features

## Usage

``` r
cuda_cell_pca(
  counts,
  n_components = 30L,
  n_hvg = 2000L,
  scale. = TRUE,
  device = c("auto", "cuda", "cpu")
)
```

## Arguments

- counts:

  Feature-by-cell matrix.

- n_components:

  Number of principal components.

- n_hvg:

  Number of highly variable features.

- scale.:

  Whether to scale selected features before PCA.

- device:

  Device passed to
  [`cudalearnr::cuda_pca()`](https://cudaverse.github.io/cudalearnr/reference/cuda_pca.html).

## Value

A `cuda_pca` object with an additional `features` element. Cell names
label PCA score rows and selected feature names label loadings.

## Examples

``` r
set.seed(2)
counts <- matrix(rpois(30 * 20, lambda = 2), 30, 20)
cuda_cell_pca(counts, n_components = 3, n_hvg = 10, device = "cpu")
#> <cuda_pca components=3 device=cpu>
#>                    PC1          PC2         PC3
#> feature_6   0.02608630  0.512558175  0.19860552
#> feature_7  -0.36704292  0.393908021  0.27650371
#> feature_13  0.40757497 -0.209287751 -0.22659441
#> feature_12  0.01657213 -0.420012007  0.54493624
#> feature_10 -0.27037240  0.199740313  0.07087489
#> feature_9  -0.19268604 -0.054987667 -0.67505649
#> feature_8  -0.34508242 -0.320837803  0.05396706
#> feature_17  0.23051134  0.005367556  0.14126006
#> feature_21 -0.23326910 -0.457835356  0.16700949
#> feature_20 -0.60112065 -0.080064874 -0.15594065
```

# Normalize a single-cell count matrix

Counts are expected with features in rows and cells in columns.

## Usage

``` r
cuda_normalize_counts(counts, scale_factor = 10000, log1p = TRUE)
```

## Arguments

- counts:

  A base matrix, `Matrix` sparse matrix, or `cudasparse`.

- scale_factor:

  Target library size.

- log1p:

  Whether to apply [`log1p()`](https://rdrr.io/r/base/Log.html) to
  non-zero normalized values.

## Value

A sparse `Matrix::dgCMatrix` with the input feature and cell names.

## Examples

``` r
counts <- Matrix::Matrix(matrix(c(1, 0, 3, 2, 4, 1), 3), sparse = TRUE)
cuda_normalize_counts(counts)
#> 3 x 2 sparse Matrix of class "dgCMatrix"
#>                       
#> [1,] 7.824446 7.957927
#> [2,] .        8.650900
#> [3,] 8.922792 7.265130
```

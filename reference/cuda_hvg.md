# Select highly variable features

Select highly variable features

## Usage

``` r
cuda_hvg(counts, n_top = 2000L, min_mean = 0)
```

## Arguments

- counts:

  Feature-by-cell count or normalized expression matrix.

- n_top:

  Number of features to select.

- min_mean:

  Minimum feature mean.

## Value

A data frame ordered by dispersion with feature statistics and a
`selected` indicator. The `index` column retains the original feature
position even when feature names are duplicated.

## Examples

``` r
set.seed(1)
counts <- Matrix::rsparsematrix(20, 10, density = 0.3)
counts@x <- abs(round(counts@x * 5))
cuda_hvg(counts, n_top = 5)
#>       feature mean   variance dispersion selected index
#> 17 feature_17  1.1 12.1000000 11.0000000     TRUE    17
#> 6   feature_6  2.0 14.4444444  7.2222222     TRUE     6
#> 12 feature_12  0.7  4.9000000  7.0000000     TRUE    12
#> 3   feature_3  1.2  7.2888889  6.0740741     TRUE     3
#> 9   feature_9  0.7  3.5666667  5.0952381     TRUE     9
#> 16 feature_16  0.7  3.5666667  5.0952381    FALSE    16
#> 1   feature_1  1.0  4.8888889  4.8888889    FALSE     1
#> 8   feature_8  0.9  3.6555556  4.0617284    FALSE     8
#> 10 feature_10  0.9  3.6555556  4.0617284    FALSE    10
#> 14 feature_14  1.4  5.6000000  4.0000000    FALSE    14
#> 4   feature_4  2.5  9.8333333  3.9333333    FALSE     4
#> 15 feature_15  0.5  1.6111111  3.2222222    FALSE    15
#> 2   feature_2  1.9  5.6555556  2.9766082    FALSE     2
#> 5   feature_5  1.7  4.9000000  2.8823529    FALSE     5
#> 7   feature_7  1.3  3.5666667  2.7435897    FALSE     7
#> 13 feature_13  0.9  2.3222222  2.5802469    FALSE    13
#> 11 feature_11  0.4  0.9333333  2.3333333    FALSE    11
#> 19 feature_19  0.4  0.4888889  1.2222222    FALSE    19
#> 20 feature_20  0.2  0.1777778  0.8888889    FALSE    20
#> 18 feature_18  0.0  0.0000000  0.0000000    FALSE    18
```

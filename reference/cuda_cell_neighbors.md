# Construct a nearest-neighbour representation

Construct a nearest-neighbour representation

## Usage

``` r
cuda_cell_neighbors(
  embedding,
  k = 15L,
  metric = c("euclidean", "cosine"),
  device = c("auto", "cuda", "cpu"),
  batch_size = 256L
)
```

## Arguments

- embedding:

  Cell-by-component matrix or `cuda_pca` result.

- k:

  Number of neighbours.

- metric:

  Distance metric.

- device:

  Computation device.

- batch_size:

  Maximum query rows per exact kNN distance block. Smaller values reduce
  peak distance memory.

## Value

A `cuda_knn` object retaining cell names from the embedding.

## Examples

``` r
cuda_cell_neighbors(
  matrix(rnorm(60), 20, 3),
  k = 5,
  batch_size = 4,
  device = "cpu"
)
#> <cuda_knn observations=20 k=5 metric=euclidean distance_device=cpu compute=cpu backend=base>
```

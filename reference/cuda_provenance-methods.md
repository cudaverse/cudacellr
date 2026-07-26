# Inspect provenance stored on native single-cell containers

These methods expose the shared
[`cudatensr::cuda_provenance()`](https://cudaverse.github.io/cudatensr/reference/cuda_provenance.html)
contract for `SingleCellExperiment` and `Seurat` objects produced by
`cudacellr`.

## Usage

``` r
# S3 method for class 'SingleCellExperiment'
cuda_provenance(x)

# S3 method for class 'Seurat'
cuda_provenance(x)
```

## Arguments

- x:

  A supported native single-cell container.

## Value

A `cuda_provenance` data frame.

# Inspect actual compute provenance

This is the shared
[`cudatensr::cuda_provenance()`](https://cudaverse.github.io/cudatensr/reference/cuda_provenance.html)
inspector, re-exposed for single-cell results.

## Usage

``` r
cuda_provenance(x)
```

## Arguments

- x:

  A cudaverse result or named list of compute stages.

## Value

A `cuda_provenance` data frame.

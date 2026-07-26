# Backend provenance and CUDA diagnostics

`cudacellr` combines sparse CPU preprocessing with CUDA-aware PCA and
nearest neighbours. A function accepting `device = "cuda"` is therefore
not necessarily an end-to-end GPU workflow. Version 0.2.0 records each
actual stage so hybrid execution remains visible.

The main tutorial below is a complete CPU workflow. The later CUDA
section is optional and executes only when diagnostics prove that a
usable device exists.

## A reproducible CPU workflow

`cudacellr` accepts features in rows and cells in columns. Start with a
small sparse count matrix:

``` r

library(cudacellr)
library(Matrix)

set.seed(42)
counts <- matrix(
  rpois(60 * 30, lambda = 2),
  nrow = 60,
  ncol = 30,
  dimnames = list(
    paste0("gene_", seq_len(60)),
    paste0("cell_", seq_len(30))
  )
)
counts[1, ] <- counts[1, ] + 1
counts <- Matrix(counts, sparse = TRUE)
```

Every stage is independently callable:

``` r

normalized <- cuda_normalize_counts(counts)
variable <- cuda_hvg(normalized, n_top = 20)
pca <- cuda_cell_pca(
  counts,
  n_hvg = 20,
  n_components = 5,
  device = "cpu"
)
neighbors <- cuda_cell_neighbors(
  pca,
  k = 5,
  batch_size = 8,
  device = "cpu"
)
```

The convenience workflow produces the same sequence with one call:

``` r

fit <- cudacell_workflow(
  counts,
  n_hvg = 20,
  n_components = 5,
  k = 5,
  batch_size = 8,
  device = "cpu"
)

fit
#> <cudacell_workflow features=60 cells=30 hvg=20 components=5 k=5 pca_device=cpu compute=cpu>
```

Normalization remains sparse, while PCA scores and neighbour results
retain cell identifiers:

``` r

class(fit$normalized)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
identical(rownames(fit$pca$x), colnames(counts))
#> [1] TRUE
identical(rownames(fit$neighbors$index), colnames(counts))
#> [1] TRUE
```

## Read the provenance table

[`cuda_provenance()`](https://cudaverse.github.io/cudatensr/reference/cuda_provenance.html)
returns the ordered stages that produced an object. The workflow table
includes fixed CPU preprocessing and the delegated `cudalearnr` PCA and
kNN stages:

``` r

cuda_provenance(normalized)
#> <cuda_provenance schema=cudaverse-stage/1 stages=1 compute=cpu>
#>          stage requested_device device backend   selection_reason fallback
#>  normalization        fixed-cpu    cpu  Matrix algorithm_cpu_only    FALSE
#>  output_device
#>            cpu
cuda_provenance(variable)
#> <cuda_provenance schema=cudaverse-stage/1 stages=2 compute=cpu>
#>          stage requested_device device backend   selection_reason fallback
#>  normalization        fixed-cpu    cpu  Matrix algorithm_cpu_only    FALSE
#>            hvg        fixed-cpu    cpu  Matrix algorithm_cpu_only    FALSE
#>  output_device
#>            cpu
#>            cpu
cuda_provenance(pca)
#> <cuda_provenance schema=cudaverse-stage/1 stages=4 compute=cpu>
#>              stage requested_device device backend   selection_reason fallback
#>      normalization        fixed-cpu    cpu  Matrix algorithm_cpu_only    FALSE
#>                hvg        fixed-cpu    cpu  Matrix algorithm_cpu_only    FALSE
#>  pca_preprocessing              cpu    cpu   stats       explicit_cpu    FALSE
#>  pca_decomposition              cpu    cpu   stats       explicit_cpu    FALSE
#>  output_device
#>            cpu
#>            cpu
#>            cpu
#>            cpu
cuda_provenance(neighbors)
#> <cuda_provenance schema=cudaverse-stage/1 stages=6 compute=cpu>
#>                   stage requested_device device backend   selection_reason
#>           normalization        fixed-cpu    cpu  Matrix algorithm_cpu_only
#>                     hvg        fixed-cpu    cpu  Matrix algorithm_cpu_only
#>       pca_preprocessing              cpu    cpu   stats       explicit_cpu
#>       pca_decomposition              cpu    cpu   stats       explicit_cpu
#>            knn_distance              cpu    cpu    base       explicit_cpu
#>  knn_neighbor_selection        fixed-cpu    cpu    base algorithm_cpu_only
#>  fallback output_device
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
workflow_provenance <- cuda_provenance(fit)
workflow_provenance
#> <cuda_provenance schema=cudaverse-stage/1 stages=6 compute=cpu>
#>                   stage requested_device device backend   selection_reason
#>           normalization        fixed-cpu    cpu  Matrix algorithm_cpu_only
#>                     hvg        fixed-cpu    cpu  Matrix algorithm_cpu_only
#>       pca_preprocessing              cpu    cpu   stats       explicit_cpu
#>       pca_decomposition              cpu    cpu   stats       explicit_cpu
#>            knn_distance              cpu    cpu    base       explicit_cpu
#>  knn_neighbor_selection        fixed-cpu    cpu    base algorithm_cpu_only
#>  fallback output_device
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
attr(workflow_provenance, "compute_device")
#> [1] "cpu"
```

The columns answer different questions:

| Column | Meaning |
|----|----|
| `requested_device` | What selected the stage: an explicit or automatic request, a fixed CPU stage, or an inherited choice. |
| `device` | Where that stage actually computed: `"cpu"` or `"cuda"`. |
| `backend` | The concrete implementation, such as `Matrix`, `base`, `stats`, or `torch`. |
| `selection_reason` | Why the stage used that device, including fixed CPU algorithms, input transfers, or CUDA availability. |
| `fallback` | `TRUE` only when an `"auto"` request selected CPU because CUDA was unavailable. |
| `output_device` | Where the stage result was materialized; public single-cell results are ordinary CPU-side R objects. |

The aggregate `compute_device` is `"cpu"` when all stages are CPU,
`"cuda"` when all stages compute on CUDA, and `"hybrid"` when a workflow
contains both. The returned object’s CPU location does not imply
fallback: a CUDA stage can intentionally materialize its result for R
and record `output_device = "cpu"`.

## Stage-by-stage backend semantics

| Function | Stage | Device/backend semantics |
|----|----|----|
| [`cuda_normalize_counts()`](https://cudaverse.github.io/cudacellr/reference/cuda_normalize_counts.md) | `normalization` | Fixed CPU `Matrix` sparse operations. |
| [`cuda_hvg()`](https://cudaverse.github.io/cudacellr/reference/cuda_hvg.md) | `hvg` | Fixed CPU `Matrix` feature statistics and base-R ranking. |
| [`cuda_cell_pca()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_pca.md) | `normalization`, `hvg` | Fixed CPU sparse preprocessing. |
| [`cuda_cell_pca()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_pca.md) | `pca_preprocessing`, `pca_decomposition` | `stats` on CPU or `torch` on CUDA; the final R result is on CPU. |
| [`cuda_cell_neighbors()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_neighbors.md) | inherited source stages | Preserves PCA or other input provenance. |
| [`cuda_cell_neighbors()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_neighbors.md) | `knn_distance` | `base` on CPU or `torch` on CUDA. |
| [`cuda_cell_neighbors()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_neighbors.md) | `knn_neighbor_selection` | Fixed CPU deterministic ordering and selection. |
| [`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md) | all stages above | CPU for sparse preprocessing; hybrid when PCA/kNN distance uses CUDA. |

An upstream `cudasparse` input contributes prefixed source stages. If a
CUDA sparse input must become a CPU `dgCMatrix`, provenance also records
an `input_materialization` transfer rather than hiding it.

## Runtime diagnostics and device requests

The shared diagnostic does not install or download torch:

``` r

diagnostics <- cudatensr::cuda_diagnostics()
diagnostics
#> <cuda_diagnostics available=FALSE devices=0 torch=not installed reason=torch_not_installed>
```

The request contract is:

- `"cpu"` guarantees the portable CPU implementation.
- `"cuda"` is strict. An unavailable runtime raises
  `cudaverse_cuda_unavailable` before CPU preprocessing starts.
- `"auto"` selects CUDA when usable; otherwise affected stages use CPU
  and record `fallback = TRUE` with a stable reason.
- Fixed CPU stages such as normalization have
  `requested_device = "fixed-cpu"` and `fallback = FALSE`. They are part
  of the algorithm, not a failed GPU request.
- A CUDA execution error after selection is reported; the workflow does
  not retry silently on CPU.

## Optional CUDA workflow

The same input can exercise every CUDA-aware single-cell stage:

``` r

if (isTRUE(diagnostics$cuda_available)) {
  cuda_fit <- cudacell_workflow(
    counts,
    n_hvg = 20,
    n_components = 5,
    k = 5,
    batch_size = 8,
    device = "cuda"
  )

  cuda_table <- cuda_provenance(cuda_fit)
  cuda_table
  attr(cuda_table, "compute_device")
} else {
  message("CUDA example skipped: ", diagnostics$reason)
}
#> CUDA example skipped: torch_not_installed
```

A successful explicit CUDA workflow is still `hybrid`: normalization,
HVG, neighbour selection, and R result materialization remain on CPU,
while PCA and kNN distance stages run on CUDA.

## Dense PCA and kNN memory

Normalization and HVG calculations retain sparse `Matrix` storage.
[`cuda_cell_pca()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_pca.md)
then selects `n_hvg` rows and materializes a dense cell-by-feature
matrix before PCA. A double copy of that selected matrix alone uses
roughly `8 * cells * n_hvg` bytes; centring, scaling, decomposition, and
CPU/GPU copies need additional working memory. Lower `n_hvg` before
lowering `n_components` when dense input memory is the constraint,
because the selected matrix is allocated before decomposition.

Exact kNN remains quadratic in time, but it does not retain a complete
cell-by-cell distance matrix. The resident double distance block is
approximately `8 * min(batch_size, cells) * cells` bytes, and the
returned distances use about `8 * cells * k` bytes plus an integer index
matrix. Lower `batch_size` to reduce peak memory; doing so does not
change the neighbours.

## Hardware-enforced parity gate

Normal R package checks must work without CUDA and may conditionally
skip the optional example above. The cudaverse hardware contract is
separate: it runs on a self-hosted runner labelled `cuda`, sets
`CUDAVERSE_REQUIRE_CUDA=true`, verifies `nvidia-smi`, requires a working
torch CUDA runtime and at least one device, and compares CPU with CUDA
PCA, kNN, and the complete single-cell workflow. It also asserts that
CUDA stages are present in provenance and the aggregate workflow is
hybrid.

The same marker can enforce a local or downstream hardware job:

``` r

require_cuda <- identical(
  tolower(Sys.getenv("CUDAVERSE_REQUIRE_CUDA", unset = "false")),
  "true"
)
if (require_cuda && !isTRUE(diagnostics$cuda_available)) {
  stop("CUDA hardware is required, but no usable device was detected.")
}
```

The package CUDA workflow runs manually, or automatically when the
repository variable `CUDAVERSE_NVIDIA_CI` is `enabled`. A required
hardware job with missing CUDA fails; it cannot pass by skipping GPU
assertions.

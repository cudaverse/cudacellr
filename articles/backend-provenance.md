# Backend provenance and CUDA diagnostics

Use `cudacellr` to run single-cell PCA and exact nearest-neighbour
search on an NVIDIA GPU, then inspect the results through ordinary R
objects. This guide shows how to prepare CUDA, run the workflow, and
check which stages used the GPU. Normalization and highly variable
feature (HVG) selection remain sparse CPU operations, so the complete
workflow is deliberately hybrid.

## Prepare CUDA

The native backend is included in `cudaverse`; you do not need torch or
LibTorch. A Windows or Linux GPU machine needs an NVIDIA driver, cuBLAS
12, and cuSOLVER 11. Follow the [platform setup
guide](https://cudaverse.github.io/cudaverse/articles/gpu-setup.html) to
install and verify these libraries. Current CUDA execution is
unavailable on macOS.

``` r

# install.packages("pak")
pak::pak("cudaverse/cudacellr")
library(cudacellr)

health <- cudaverse::cuda_diagnostics()
health$summary
health$next_steps
health$selected_backend  # "native" for the lightweight CUDA path
cudaverse::cuda_select_device("cuda")
```

Run the CUDA chunks below on your configured GPU machine. They are
displayed without execution when this article is built so installing or
checking the R package does not require a GPU. The small CPU reference
at the end is executed during the article build.

## Prepare a count matrix

`cudacellr` expects features in rows and cells in columns. This small
named matrix makes the example reproducible; replace it with your own
count matrix for analysis.

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

## Normalize, select features, and run PCA and kNN

One call reuses each preprocessing stage and returns normalized
expression, feature statistics, the PCA model, and exact neighbours:

``` r

fit <- cudacell_workflow(
  counts,
  n_hvg = 20,
  n_components = 5,
  k = 5,
  batch_size = 8,
  device = "cuda"
)

fit
head(fit$pca$x)
head(fit$neighbors$index)
head(fit$neighbors$distance)

identical(rownames(fit$pca$x), colnames(counts))
identical(rownames(fit$neighbors$index), colnames(counts))
```

You can also call stages independently. For example, pass the PCA model
directly to neighbour search:

``` r

pca <- cuda_cell_pca(
  counts,
  n_hvg = 20,
  n_components = 5,
  device = "cuda"
)
neighbors <- cuda_cell_neighbors(
  pca,
  k = 5,
  batch_size = 8,
  device = "cuda"
)
```

Native PCA returns ordinary R score and loading matrices while retaining
a device-side score cache. Passing the unmodified PCA result directly to
kNN allows reuse of that cache without uploading the scores again.
Distance blocks and deterministic top-k selection run on CUDA; the final
neighbour indices and distances are returned as R matrices. The public
PCA model is also materialized on the host, so this is not a promise
that only the final neighbour output ever crosses the device boundary.

## Read the provenance table

[`cuda_provenance()`](https://cudaverse.github.io/cudaverse/reference/cuda_provenance.html)
returns the ordered compute stages that produced an object, including
inherited preprocessing:

``` r

workflow_provenance <- cuda_provenance(fit)
workflow_provenance
attr(workflow_provenance, "compute_device")
```

| Column | Meaning |
|----|----|
| `requested_device` | An explicit or automatic request, a fixed CPU stage, or an inherited choice. |
| `device` | Where that stage computed: `"cpu"` or `"cuda"`. |
| `backend` | The implementation, such as `Matrix`, `base`, `stats`, `native`, or the optional `torch` backend. |
| `selection_reason` | Why that device was used, including fixed CPU algorithms, input transfers, or CUDA availability. |
| `fallback` | Whether an `"auto"` request used CPU because CUDA was unavailable. |
| `output_device` | Where the stage result was materialized; a CUDA calculation can return an R object on CPU. |

The aggregate `compute_device` is `"cpu"` when all recorded stages use
CPU, `"cuda"` when all use CUDA, and `"hybrid"` when they contain both.
A successful `cudacell_workflow(..., device = "cuda")` is hybrid because
normalization and HVG selection run on CPU. Returning a result to R does
not itself mean the calculation fell back to CPU.

## Native stage boundaries

These semantics describe the native backend in `cudaverse` 0.4.1:

| Function | Stage | Device/backend semantics |
|----|----|----|
| [`cuda_normalize_counts()`](https://cudaverse.github.io/cudacellr/reference/cuda_normalize_counts.md) | `normalization` | CPU `Matrix` sparse operations. |
| [`cuda_hvg()`](https://cudaverse.github.io/cudacellr/reference/cuda_hvg.md) | `hvg` | CPU `Matrix` feature statistics and base-R ranking. |
| [`cuda_cell_pca()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_pca.md) | `normalization`, `hvg` | CPU sparse preprocessing; selected features are materialized as a dense input matrix. |
| [`cuda_cell_pca()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_pca.md) | `pca_preprocessing`, `pca_decomposition` | CUDA centring/scaling and decomposition with `native`; the public PCA model is returned to R. |
| [`cuda_cell_pca()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_pca.md) | `pca_scores_resident` | A retained CUDA score cache supports subsequent native operations. |
| [`cuda_cell_neighbors()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_neighbors.md) | inherited source stages | Preserves PCA or other input provenance. |
| [`cuda_cell_neighbors()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_neighbors.md) | `knn_distance`, `knn_neighbor_selection` | CUDA distance blocks and deterministic top-k with `native`; final neighbour matrices are returned to R. |
| [`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md) | all stages above | Hybrid CPU preprocessing and CUDA PCA/kNN. |

The optional torch compatibility backend can have different boundaries:
device-side neighbour selection depends on its stable-sort support, and
older torch APIs may use a recorded CPU selection path. Neither torch
nor LibTorch is needed for the native workflow shown here. Inspect the
actual provenance table rather than assuming that every CUDA backend is
identical.

An upstream `cudasparse` input contributes prefixed source stages.
Single-cell normalization still needs a CPU `dgCMatrix`; materializing a
CUDA sparse input therefore records an `input_materialization` transfer.

## Device requests and failures

- `"cuda"` is strict: an unavailable runtime raises
  `cudaverse_cuda_unavailable` before single-cell preprocessing starts.
- `"auto"` uses CUDA when eligible; otherwise the requested numerical
  stages run on CPU and record fallback with a reason.
- `"cpu"` selects CPU computation for a reference run or a machine
  without CUDA.
- Fixed CPU stages such as normalization have
  `requested_device = "fixed-cpu"` and `fallback = FALSE`. They are part
  of the workflow, not a failed GPU request.
- A CUDA execution error after selection is reported; the workflow does
  not retry silently on CPU.

After changing drivers or CUDA library paths, restart R and run
[`cudaverse::cuda_diagnostics()`](https://cudaverse.github.io/cudaverse/reference/cuda_diagnostics.html)
again. Its `summary` and `next_steps` fields help distinguish missing
runtime libraries from a failed device self-test.

## Plan memory and workload size

Normalization and HVG calculations retain sparse `Matrix` storage.
[`cuda_cell_pca()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_pca.md)
then selects `n_hvg` rows and materializes a dense cell-by-feature
matrix. One double-precision copy of this selected matrix uses roughly
`8 * cells * n_hvg` bytes; centring, scaling, decomposition, and
host/device copies require additional memory. Lower `n_hvg` when dense
input memory is the constraint: the input is allocated before
decomposition, so reducing `n_components` alone does not make that
matrix smaller.

The PCA model retains host scores/loadings and, for native CUDA, cached
device scores. Remove unused models and call
[`gc()`](https://rdrr.io/r/base/gc.html) to allow their storage to be
released. `cudaverse::cuda_memory_info("cuda")` can help inspect device
memory; whole-device usage can include other applications.

Exact kNN remains quadratic in time but uses bounded distance blocks
rather than retaining the full cell-by-cell matrix. The double distance
block uses approximately `8 * min(batch_size, cells) * cells` bytes,
plus backend workspace. Returned distances use about `8 * cells * k`
bytes plus an integer index matrix. Lower `batch_size` to reduce peak
memory; it does not change the exact neighbour definition or its
original-row-index tie breaking.

This tiny tutorial demonstrates the API, not a GPU speed advantage.
Small inputs may be dominated by startup and transfer costs. Benchmark
your complete workflow, including CPU normalization/HVG and transfers,
before making performance claims about single-cell analysis.

## A small CPU reference

The same API can produce a reference result without CUDA:

``` r

reference <- cudacell_workflow(
  counts,
  n_hvg = 20,
  n_components = 5,
  k = 5,
  batch_size = 8,
  device = "cpu"
)

reference
#> <cudacell_workflow features=60 cells=30 hvg=20 components=5 k=5 pca_device=cpu compute=cpu>
identical(rownames(reference$pca$x), colnames(counts))
#> [1] TRUE
identical(rownames(reference$neighbors$index), colnames(counts))
#> [1] TRUE
cuda_provenance(reference)
#> <cuda_provenance schema=cudaverse-stage/1 stages=6 compute=cpu>
#>                   stage requested_device device backend   selection_reason
#>           normalization        fixed-cpu    cpu  Matrix algorithm_cpu_only
#>                     hvg        fixed-cpu    cpu  Matrix algorithm_cpu_only
#>       pca_preprocessing              cpu    cpu   stats       explicit_cpu
#>       pca_decomposition              cpu    cpu   stats       explicit_cpu
#>            knn_distance              cpu    cpu    base       explicit_cpu
#>  knn_neighbor_selection              cpu    cpu    base       explicit_cpu
#>  fallback output_device
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
#>     FALSE           cpu
```

For a GPU comparison, rerun the CUDA workflow above on the same input
and compare numerical results within an appropriate tolerance. PCA
vectors may differ in sign, so compare reconstruction or subspaces
instead of raw loading signs. Check neighbour indices and distances,
including distance ties, and retain the runtime diagnostics and
provenance with the results.

Package installation and ordinary checks do not require CUDA. Hardware
validation must actually execute the CUDA calls, verify the selected
backend and recorded stages, and check numerical agreement; a skipped
GPU example or job is not evidence of a successful hardware test.

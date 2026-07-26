# A native SingleCellExperiment workflow

[`cudacell_sce()`](https://cudaverse.github.io/cudacellr/reference/cudacell_sce.md)
runs the complete cudacellr preprocessing, PCA, and nearest neighbour
workflow while keeping results in the native `SingleCellExperiment` data
model. It returns a modified copy; the input object is not changed.

The Bioconductor packages are optional so matrix-only users do not
acquire a large mandatory dependency tree:

``` r

BiocManager::install("SingleCellExperiment")
```

## Start with a populated object

This example intentionally includes an unrelated assay, row and column
metadata, a pre-existing reduced dimension, top-level metadata, size
factors, and an alternative experiment. All of them must survive the
workflow.

``` r

library(cudacellr)
library(Matrix)
library(SingleCellExperiment)
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:Matrix':
#> 
#>     expand, unname
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians

set.seed(42)
counts <- matrix(
  rpois(80 * 30, lambda = 2),
  nrow = 80,
  ncol = 30,
  dimnames = list(
    paste0("gene_", seq_len(80)),
    paste0("cell_", seq_len(30))
  )
)
counts[1, ] <- counts[1, ] + 1
counts <- Matrix(counts, sparse = TRUE)

sce <- SingleCellExperiment(
  assays = list(
    counts = counts,
    untouched = counts * 2
  ),
  rowData = S4Vectors::DataFrame(
    symbol = paste0("symbol_", seq_len(nrow(counts)))
  ),
  colData = S4Vectors::DataFrame(
    batch = rep(c("one", "two"), each = 15)
  ),
  metadata = list(owner = "unchanged"),
  reducedDims = list(
    EXISTING = matrix(rnorm(60), nrow = 30)
  )
)
sizeFactors(sce) <- rep(1, ncol(sce))
altExp(sce, "spike") <- SingleCellExperiment(
  assays = list(
    counts = matrix(
      rpois(5 * ncol(sce), 1),
      nrow = 5,
      dimnames = list(paste0("spike_", 1:5), colnames(sce))
    )
  )
)
```

## Run the workflow

The default output names are namespaced so common pre-existing
`logcounts` and `PCA` entries are not overwritten:

``` r

result <- cudacell_sce(
  sce,
  assay = "counts",
  n_hvg = 30,
  n_components = 5,
  k = 5,
  batch_size = 8,
  device = "cpu"
)

assayNames(result)
#> [1] "counts"             "untouched"          "cudacell_logcounts"
reducedDimNames(result)
#> [1] "EXISTING"     "CUDACELL_PCA"
colPairNames(result)
#> [1] "CUDACELL_KNN"
```

The result locations follow their SingleCellExperiment semantics:

``` r

assay(result, "cudacell_logcounts")[1:3, 1:3]
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>          cell_1   cell_2   cell_3
#> gene_1 5.652808 5.208492 5.773037
#> gene_2 5.430541 4.120760 4.861401
#> gene_3 4.057303 4.120760 4.861401
reducedDim(result, "CUDACELL_PCA")[1:3, , drop = FALSE]
#>              PC1      PC2       PC3        PC4         PC5
#> cell_1 -1.010505 2.472902 -1.934263  1.0690482 -0.05279964
#> cell_2 -1.130595 1.670652  1.065843 -2.2639677 -3.15504850
#> cell_3  2.060318 3.225506 -1.481663  0.6235578 -0.62236913
rowData(result)[
  1:3,
  c(
    "cudacell_mean",
    "cudacell_variance",
    "cudacell_dispersion",
    "cudacell_hvg",
    "cudacell_hvg_rank"
  )
]
#> DataFrame with 3 rows and 5 columns
#>        cudacell_mean cudacell_variance cudacell_dispersion cudacell_hvg
#>            <numeric>         <numeric>           <numeric>    <logical>
#> gene_1       5.24144          0.289807           0.0552916        FALSE
#> gene_2       4.51740          1.797625           0.3979330        FALSE
#> gene_3       4.02805          2.814162           0.6986410        FALSE
#>        cudacell_hvg_rank
#>                <integer>
#> gene_1                80
#> gene_2                73
#> gene_3                50
colPair(result, "CUDACELL_KNN")
#> SelfHits object with 150 hits and 2 metadata columns:
#>              from        to |  distance      rank
#>         <integer> <integer> | <numeric> <integer>
#>     [1]         1         8 |   2.57148         4
#>     [2]         1         9 |   1.02935         1
#>     [3]         1        17 |   2.44363         3
#>     [4]         1        18 |   1.96267         2
#>     [5]         1        29 |   2.93467         5
#>     ...       ...       ... .       ...       ...
#>   [146]        30         7 |   1.87794         1
#>   [147]        30        19 |   3.20619         5
#>   [148]        30        20 |   2.29256         2
#>   [149]        30        24 |   2.95473         3
#>   [150]        30        27 |   3.15710         4
#>   -------
#>   nnode: 30
```

`CUDACELL_KNN` is a directed `SelfHits` object. Its metadata columns
contain distance and neighbour rank. Unlike an integer matrix hidden in
top-level metadata, this representation automatically remaps cell
indices when the SCE is subset.

## Verify preservation

The original object is unchanged and unrelated fields are retained:

``` r

identical(assay(sce, "counts"), assay(result, "counts"))
#> [1] TRUE
identical(assay(sce, "untouched"), assay(result, "untouched"))
#> [1] TRUE
identical(rowData(sce)$symbol, rowData(result)$symbol)
#> [1] TRUE
identical(colData(sce)$batch, colData(result)$batch)
#> [1] TRUE
identical(reducedDim(sce, "EXISTING"), reducedDim(result, "EXISTING"))
#> [1] TRUE
identical(altExp(sce, "spike"), altExp(result, "spike"))
#> [1] TRUE
identical(sizeFactors(sce), sizeFactors(result))
#> [1] TRUE
identical(metadata(sce)$owner, metadata(result)$owner)
#> [1] TRUE
```

Feature statistics are mapped by their original integer row position,
not by feature name. Duplicate feature names are therefore safe.

## Normalization and size factors

The default normalized assay uses natural `log1p`, so it is named
`cudacell_logcounts` rather than the Bioconductor convention
`logcounts`, which commonly implies log base 2. Set `log1p = FALSE` to
create `cudacell_normalized`, or choose an explicit output name:

``` r

result <- cudacell_sce(
  sce,
  log1p = FALSE,
  normalized_assay = "normalized_counts",
  scale_factor = 10000,
  n_hvg = 30,
  n_components = 5,
  k = 5,
  device = "cpu"
)
```

The actual divisor, `library_size / scale_factor`, is always stored in
`colData(result)$cudacell_size_factor`. Canonical `sizeFactors(result)`
are left untouched unless `set_size_factors = TRUE`. Replacing existing
canonical size factors additionally requires `overwrite = TRUE`.

## Collision and memory safety

Every target assay, reduced dimension, colPair, rowData column, colData
column, and the `cudacellr` metadata record is checked before
computation. Rerunning with the same names fails rather than discarding
data:

``` r

cudacell_sce(
  result,
  n_hvg = 30,
  n_components = 5,
  k = 5,
  device = "cpu"
)
#> Error:
#> ! cudacellr output fields already exist: assay 'cudacell_logcounts', reducedDim 'CUDACELL_PCA', colPair 'CUDACELL_KNN', rowData 'cudacell_mean', rowData 'cudacell_variance', rowData 'cudacell_dispersion', rowData 'cudacell_hvg', rowData 'cudacell_hvg_rank', colData 'cudacell_size_factor', metadata 'cudacellr'. Set `overwrite = TRUE` to replace only these fields.
```

Use `overwrite = TRUE` only when replacing those exact cudacellr outputs
is intentional. No unrelated field is replaced.

On-disk and other delayed assays are also rejected by default. cudacellr
does not silently turn a potentially enormous delayed assay into an
in-memory matrix. `realize = TRUE` is an explicit opt-in after checking
the expected memory requirement.

## Provenance and downstream embeddings

The lightweight metadata record contains parameters, output locations,
and the same validated `cudaverse-stage/1` compute stages as the matrix
workflow:

``` r

cuda_provenance(result)
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
metadata(result)$cudacellr$parameters
#> $n_hvg
#> [1] 30
#> 
#> $n_components
#> [1] 5
#> 
#> $k
#> [1] 5
#> 
#> $requested_device
#> [1] "cpu"
#> 
#> $batch_size
#> [1] 8
#> 
#> $scale_factor
#> [1] 10000
#> 
#> $log1p
#> [1] TRUE
#> 
#> $min_mean
#> [1] 0
#> 
#> $scale
#> [1] TRUE
```

Per-feature, per-cell, and graph-shaped results do not live in top-level
metadata because that metadata is not remapped by row or column
subsetting. For a downstream embedding, pass `result` directly to a
compatible `cudaembedr` function or explicitly use
`reducedDim(result, "CUDACELL_PCA")`.

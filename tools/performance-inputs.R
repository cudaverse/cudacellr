# Reproducible inputs; run outside R CMD check. No download occurs here.
args <- commandArgs(TRUE)
root <- normalizePath(args[[1]], mustWork = TRUE)
library(Matrix)
library(jsonlite)
library(digest)
archive <- file.path(root, "pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz")
expected <- "060af9e4370c6ad877e2987cf8cdc63a6d5c84e57875580f2be6c8e90748cf63"
stopifnot(identical(digest(archive, algo = "sha256", file = TRUE), expected))
files <- untar(archive, list = TRUE)
stopifnot(!any(grepl("(^/|(^|/)\\.\\.(/|$))", files)))
untar(archive, exdir = root)
cases <- list()
emit <- function(id, x, hvg, components, seed = NA_integer_) {
  x <- as(x, "dgCMatrix")
  path <- file.path(root, paste0(id, ".rds"))
  saveRDS(x, path, compress = FALSE)
  cases[[length(cases) + 1L]] <<- list(
    id = id, file = basename(path), sha256 = digest(path, algo = "sha256", file = TRUE),
    cells = ncol(x), features = nrow(x), nnz = length(x@x), seed = seed,
    hvg = hvg, components = components, k = 15L, batch_size = 256L)
}
for (spec in list(c(1000, 50, 1), c(10000, 100, 42), c(50000, 128, 2026))) {
  set.seed(spec[3])
  x <- matrix(rpois(spec[1] * spec[2], 0.3), nrow = spec[2])
  x[1, ] <- x[1, ] + 1L
  dimnames(x) <- list(paste0("feature_", seq_len(nrow(x))),
                     paste0("cell_", seq_len(ncol(x))))
  emit(paste0("legacy_", spec[1]), Matrix(x, sparse = TRUE),
       floor(spec[2] * 0.8), 30L, spec[3])
}
# Sparse stress fixture: 10,000 genes, 50,000 cells, ~0.5% nonzeros.
set.seed(20261019)
x <- rsparsematrix(10000, 50000, density = 0.005,
                   rand.x = function(n) rpois(n, 2) + 1)
x[1, ] <- x[1, ] + 1
dimnames(x) <- list(paste0("gene_", seq_len(nrow(x))),
                   paste0("cell_", seq_len(ncol(x))))
emit("sparse_50000", x, 2000L, 50L, 20261019L)
folder <- file.path(root, "filtered_feature_bc_matrix")
x <- readMM(gzfile(file.path(folder, "matrix.mtx.gz")))
features <- read.delim(gzfile(file.path(folder, "features.tsv.gz")), header = FALSE)
barcodes <- read.delim(gzfile(file.path(folder, "barcodes.tsv.gz")), header = FALSE)
dimnames(x) <- list(features[[1]], barcodes[[1]])
keep <- if (ncol(features) >= 3L) features[[3]] == "Gene Expression" else rep(TRUE, nrow(x))
x <- x[keep, , drop = FALSE]
stopifnot(all(colSums(x) > 0))
emit("pbmc_10k_v3", x, 2000L, 50L)
write_json(list(schema = "cudacellr-performance-inputs/1",
  pbmc = list(url = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz",
              sha256 = expected, license = "CC BY 4.0", attribution = "10x Genomics"),
  cases = cases), file.path(root, "inputs.json"), auto_unbox = TRUE, pretty = TRUE)

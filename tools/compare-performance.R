args <- commandArgs(TRUE)
library(jsonlite)
root <- args[1]
cases <- c("legacy_1000", "legacy_10000", "legacy_50000", "sparse_50000", "pbmc_10k_v3")
close <- function(a, b, rtol = 1e-8, atol = 1e-10) {
  identical(dim(a), dim(b)) && all(is.finite(a)) && all(is.finite(b)) &&
    all(abs(a - b) <= atol + rtol * abs(b))
}
results <- list()
for (id in cases) {
  paths <- file.path(root, paste0(id, "-", c("baseline", "candidate"), "-timing.json"))
  a <- readRDS(paste0(paths[1], ".result.rds")); b <- readRDS(paste0(paths[2], ".result.rds"))
  x <- read_json(paths[1], simplifyVector = TRUE); y <- read_json(paths[2], simplifyVector = TRUE)
  recon <- TRUE
  for (start in seq.int(1L, nrow(a$scores), 512L)) {
    rows <- start:min(nrow(a$scores), start + 511L)
    recon <- recon && close(tcrossprod(b$scores[rows, , drop = FALSE], b$rotation),
                           tcrossprod(a$scores[rows, , drop = FALSE], a$rotation))
  }
  checks <- c(normalization = identical(a$normalized_hash, b$normalized_hash),
    hvg = identical(a$hvg_hash, b$hvg_hash), features = identical(a$features, b$features),
    identifiers = identical(dimnames(a$scores), dimnames(b$scores)) &&
      identical(dimnames(a$rotation), dimnames(b$rotation)) &&
      identical(dimnames(a$index), dimnames(b$index)),
    sdev = close(a$sdev, b$sdev), center = close(a$center, b$center), scale = close(a$scale, b$scale),
    projector = close(tcrossprod(b$rotation), tcrossprod(a$rotation), 0, 1e-8),
    reconstruction = recon, neighbors = identical(a$index, b$index),
    distances = close(b$distance, a$distance),
    same_environment = identical(x$host, y$host) && identical(x$r_version, y$r_version) &&
      identical(x$core_version, y$core_version) && identical(x$case, y$case) &&
      identical(x$input_sha256, y$input_sha256),
    protocol = x$warmups == 5 && y$warmups == 5 && x$samples == 10 && y$samples == 10,
    cleanup = x$gpu_after$allocated_bytes == x$gpu_before$allocated_bytes &&
      y$gpu_after$allocated_bytes == y$gpu_before$allocated_bytes)
  ratio <- median(y$elapsed_seconds) / median(x$elapsed_seconds)
  reduction <- 1 - y$host_peak_rss_bytes / x$host_peak_rss_bytes
  results[[id]] <- list(checks = as.list(checks), correctness_pass = all(checks),
    baseline_median_seconds = median(x$elapsed_seconds), candidate_median_seconds = median(y$elapsed_seconds),
    runtime_ratio = ratio, host_peak_reduction = reduction,
    baseline_host_peak_bytes = x$host_peak_rss_bytes, candidate_host_peak_bytes = y$host_peak_rss_bytes,
    runtime_gate = ratio <= 1.1, memory_gate = if (id == "sparse_50000") reduction >= 0.25 else NA)
}
pass <- all(vapply(results, function(x) x$correctness_pass && x$runtime_gate, logical(1))) &&
  results$sparse_50000$memory_gate
write_json(list(schema = "cudacellr-performance-comparison/1", cases = results,
  passed = pass, shared_gpu = TRUE,
  note = "Observed shared-device timings; no uncontended speedup claim. Review process observations before integration."),
  file.path(root, "comparison.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
print(results)
if (!pass) quit(status = 1L)

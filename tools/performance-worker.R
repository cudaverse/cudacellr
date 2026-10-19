args <- commandArgs(TRUE)
stopifnot(length(args) == 5L)
data <- normalizePath(args[1]); id <- args[2]; variant <- args[3]
mode <- args[4]; out <- args[5]
library(cudacellr)
library(Matrix)
library(jsonlite)
library(digest)
options(cudaverse.cuda_backends = "native")
stopifnot(as.character(packageVersion("cudaverse")) == "0.4.2")
stopifnot(as.character(packageVersion("cudacellr")) ==
            if (variant == "baseline") "0.4.0" else "0.4.0.9000")
manifest <- read_json(file.path(data, "inputs.json"), simplifyVector = FALSE)
spec <- Filter(function(x) identical(x$id, id), manifest$cases)[[1]]
path <- file.path(data, spec$file)
stopifnot(identical(digest(path, algo = "sha256", file = TRUE), spec$sha256))
counts <- readRDS(path)
stopifnot(identical(cudaverse::cuda_select_device("cuda")$backend, "native"))
sync <- getFromNamespace(".native_synchronize", "cudaverse")
rss <- function(field = "VmHWM") {
  line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", line)) * 1024
}
gpu <- function() {
  text <- system2("nvidia-smi", c("--query-compute-apps=pid,used_gpu_memory",
                                "--format=csv,noheader,nounits"), stdout = TRUE, stderr = TRUE)
  list(time = format(Sys.time(), tz = "UTC", usetz = TRUE),
       pid = Sys.getpid(), status = if (is.null(attr(text, "status"))) 0L else attr(text, "status"),
       compute_processes = text)
}
run <- function() cudacell_workflow(counts, n_hvg = spec$hvg,
  n_components = spec$components, k = spec$k, batch_size = spec$batch_size,
  scale. = TRUE, device = "cuda")
plain <- function(x) { attr(x, "cudaverse_native_state") <- NULL; x }
observations <- list(); timings <- numeric(); stages <- list()
before <- unclass(cudaverse::cuda_memory_info("cuda"))
rss_before <- rss("VmRSS")
if (mode == "timing") {
  for (i in seq_len(5L)) { invisible(run()); gc(); sync() }
  for (i in seq_len(10L)) {
    gc(); sync(); observations[[length(observations) + 1L]] <- gpu()
    start <- proc.time()[[3]]; result <- run(); sync()
    timings[i] <- proc.time()[[3]] - start
    observations[[length(observations) + 1L]] <- gpu()
    if (i == 1L) {
      saveRDS(list(normalized_hash = digest(list(result$normalized@x, result$normalized@i,
        result$normalized@p, dimnames(result$normalized)), algo = "sha256"),
        hvg_hash = digest(result$variable_features, algo = "sha256"),
        scores = plain(result$pca$x), rotation = result$pca$rotation,
        sdev = result$pca$sdev, center = result$pca$center, scale = result$pca$scale,
        features = result$pca$features, index = result$neighbors$index,
        distance = result$neighbors$distance), paste0(out, ".result.rds"), compress = FALSE)
      provenance <- as.data.frame(cudaverse::cuda_provenance(result))
    }
    result <- NULL
  }
} else if (mode == "profile") {
  measure <- function(name, expr) {
    sync(); t <- proc.time()[[3]]; value <- force(expr); sync()
    stages[[name]] <<- list(elapsed = proc.time()[[3]] - t, host_rss_bytes = rss("VmRSS"),
                           gpu = unclass(cudaverse::cuda_memory_info("cuda")))
    value
  }
  normalized <- measure("normalization", cuda_normalize_counts(counts))
  variable <- measure("feature_selection", cuda_hvg(normalized, n_top = spec$hvg))
  selected <- variable[variable$selected, , drop = FALSE]
  input <- measure("pca_preparation", {
    if (variant == "candidate") {
      getFromNamespace(".cell_pca_input", "cudacellr")(
        normalized, selected$index, selected$feature, "cuda")
    } else t(as.matrix(normalized[selected$index, , drop = FALSE]))
  })
  fit <- measure("pca_including_upload", cudaverse::cuda_pca(input,
    n_components = spec$components, scale. = TRUE, device = "cuda"))
  neighbors <- measure("neighbors", cuda_cell_neighbors(fit,
    k = spec$k, batch_size = spec$batch_size, device = "cuda"))
  provenance <- as.data.frame(cudaverse::cuda_provenance(neighbors))
  transfer <- list(input_representation = class(input)[1],
    host_pca_input_bytes = as.numeric(object.size(input)),
    payload_estimate_bytes = if (variant == "candidate")
      length(input$values) * 8 + length(input$col_index) * 4 + length(input$row_ptr) * 4 else
      length(input) * 8,
    estimate_scope = "Input numeric/CSR payload estimate; not measured driver traffic",
    scores_reused = any(grepl("scores_resident", provenance$stage)))
} else stop("Unknown mode")
peak <- unclass(cudaverse::cuda_memory_info("cuda"))
host_peak <- rss()
rm(list = intersect(c("result", "normalized", "variable", "selected", "input", "fit", "neighbors"), ls()))
gc(); sync()
after <- unclass(cudaverse::cuda_memory_info("cuda"))
write_json(list(schema = "cudacellr-performance/1", complete = TRUE, case = spec,
  variant = variant, mode = mode, r_version = R.version.string,
  package_version = as.character(packageVersion("cudacellr")), core_version = "0.4.2",
  host = Sys.info()[["nodename"]], pid = Sys.getpid(), input_sha256 = spec$sha256,
  elapsed_seconds = timings, warmups = if (mode == "timing") 5L else 0L,
  samples = length(timings), stages = stages, host_peak_rss_bytes = host_peak,
  host_rss_before_bytes = rss_before, gpu_before = before, gpu_peak = peak,
  gpu_after = after, gpu_observations = observations, provenance = provenance,
  transfers = if (exists("transfer")) transfer else NULL,
  timing_scope = "Synchronized public workflow call; GC and GPU process queries outside timing",
  host_memory_scope = "Linux process VmHWM, including R/runtime/input/results; separate fresh process per variant and mode",
  shared_gpu = TRUE), out, pretty = TRUE, auto_unbox = TRUE, na = "null")

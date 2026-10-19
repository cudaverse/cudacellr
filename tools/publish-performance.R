# Publish an explicitly labeled derivative without unrelated process identities.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, !dir.exists(args[2]))
dir.create(args[2], recursive = TRUE)
for (path in list.files(args[1], pattern = "-(timing|profile)\\.json$", full.names = TRUE)) {
  report <- jsonlite::read_json(path, simplifyVector = FALSE)
  benchmark_pid <- as.character(report$pid)
  report$pid <- NULL
  report$host <- "shared-a100-validation-host"
  report$gpu_observations <- lapply(report$gpu_observations, function(observation) {
    entries <- strsplit(unlist(observation$compute_processes), ",", fixed = TRUE)
    process_id <- vapply(entries, function(x) trimws(x[1]), character(1))
    memory <- vapply(entries, function(x) as.numeric(trimws(x[2])), numeric(1))
    list(time = observation$time, status = observation$status,
      benchmark_gpu_mib = sum(memory[process_id == benchmark_pid]),
      other_compute_process_count = sum(process_id != benchmark_pid),
      other_compute_gpu_mib = sum(memory[process_id != benchmark_pid]))
  })
  report$source_report_sha256 <- digest::digest(path, file = TRUE, algo = "sha256")
  report$publication_note <- paste(
    "Public derivative: host label and process identities removed.",
    "GPU process observations retain counts and memory totals. Original retained separately.")
  jsonlite::write_json(report, file.path(args[2], basename(path)),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
}
for (name in c("inputs.json", "comparison.json")) {
  stopifnot(file.copy(file.path(args[1], name), file.path(args[2], name)))
}

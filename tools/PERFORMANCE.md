# Sparse PCA performance protocol

Run against the unchanged submitted cudacellr 0.4.0 archive and a source-built
candidate, both using CRAN cudaverse 0.4.2. Keep separate installed libraries.
The scripts require Linux, native CUDA, Matrix, jsonlite and digest. They are
maintainer tools excluded from the CRAN archive, with no new package imports.

Download the filtered Matrix Market archive linked from the public 10x Genomics
"10k PBMCs from a Healthy Donor (v3 chemistry)" dataset. It is CC BY 4.0; retain
10x Genomics attribution. `performance-inputs.R DATA_DIRECTORY` verifies the
pinned source SHA-256, builds the three historical synthetic fixtures and a
10,000-gene/50,000-cell sparse stress fixture, and saves checksummed RDS inputs.

`run-performance.sh ROOT NEW_REPORT_DIRECTORY` runs fresh R processes for each
case, variant and mode. ROOT contains `source`, `data`, `library` (candidate and
core), and `baseline-library` (submitted companion). Five warmups and ten
synchronized measurements exercise the actual public workflow. GPU process
queries, garbage collection, serialization and profiling are outside timing.
The stage profile is a separate process: normalization, feature selection,
input preparation, PCA including upload, and neighbours. It records provenance
and clearly labeled payload estimates, not measured transfer bandwidth.

Host peak memory is Linux process VmHWM, including runtime, input and results.
GPU owned allocation peaks use cuda_memory_info; whole-device usage is also
retained and includes other clients. Existing shared GPU workloads must remain
untouched. Interpret timing as shared-device observations, never uncontended
hardware speedups. Preserve all failed attempts and report exact source,
archive, script and installed-package hashes alongside results.

`compare-performance.R REPORT_DIRECTORY` requires all fixtures, ten samples,
five warmups, exact normalization/HVG/identifiers/neighbours, existing FP64
PCA reconstruction/projector tolerances and exact native allocation cleanup.
The larger sparse fixture must reduce total process peak RSS by at least 25%.
Each case must have candidate/baseline median runtime at most 1.10. A failed
gate keeps the optimization unmerged. Review shared-GPU observations before
using the measured comparison for integration.

When only the candidate changes, an optional third driver argument names a
completed baseline report directory. The driver copies and hashes those
unchanged baseline results, then measures the new candidate. The comparison
also requires matching host, R/core versions, input hashes and parameters.
`CUDACELLR_CANDIDATE_LIBRARY` and `CUDACELLR_BENCHMARK_SOURCE` can point at
separate installations and source exports, preserving earlier attempts.

# Sparse PCA benchmark evidence — October 19, 2026

The final candidate passed all five correctness/runtime comparisons and the
large sparse memory gate. These are observations on one shared-use A100 host.
They do not establish maximum capacity or an uncontended hardware speedup.

## Identity and reproduction

- Candidate runtime: `b1b2262e453f153b6a24a5bb30bcfd031116dc17`, version
  `0.4.0.9000`. Later commits in this milestone add tooling, evidence and docs.
- Baseline: unchanged submitted cudacellr `0.4.0`, archive SHA-256
  `bb37db8a92a5b4203e5086082a6ae62ea7ae0af3397ea5c74e0872939f22ea80`.
  Original source was `02a5b1e130aaf5e8d3db4d6e65bba6eed74bfa40` plus
  three existing vignette author additions.
- Core: CRAN cudaverse `0.4.2`, published archive SHA-256
  `e31b68130556df1290403b6a3e05365dcbe0f4aa5c39e857b78d7e2c27c604c9`.
- Hardware: NVIDIA A100-SXM4-40GB, AMD EPYC 7713 VM with 32 visible CPUs.
  Ubuntu 24.04.2, R 4.3.3, Matrix 1.6-5, reference BLAS/LAPACK 3.12.0;
  `OMP_NUM_THREADS=2`, `OPENBLAS_NUM_THREADS=2`, native CUDA, float64 PCA.
- [Input manifest](inputs.json): download and derived-input hashes, seeds,
  dimensions, nonzero counts, HVG/PC/k/batch parameters.
- [Maintainer protocol](../../tools/PERFORMANCE.md): input preparation,
  isolated libraries, five warmups, ten synchronized public-workflow calls,
  stage profiles, comparison gates and separate Nsight transfer traces.

The baseline driver used tooling at `5f245a6`; final candidate timing/profile
tooling was `b1b2262`. The latter verifies matching host, R/core versions and
checksummed input parameters when reusing the unchanged baseline.
Transfer tooling was `3fb52da`, completed by `2540ad7` using the installed
matching QdstrmImporter after Debian's nsys could not locate it automatically.
No timed workflow was run under the profiler.

## Contents and interpretation

- [Final comparison](comparison.json): all numerical, cleanup, protocol,
  runtime and memory acceptance decisions. `passed` is true.
- `*-timing.json`: all ten elapsed times and process-wide peak host memory.
- `*-profile.json`: separate single-pass stage measurements, native allocation
  peaks, provenance, host input size and explicitly estimated CSR payload.
- `*-stats.csv`: measured CUDA transfer sizes/time, API and kernel summaries.
  Sizes are reported by Nsight in decimal MB. These traces include capability
  probes and return copies. The CSR estimate excludes the native COO row-index
  upload; use the trace totals for measured transfer traffic.
- [Initial comparison](first-candidate-comparison.json): preserved failed
  small-input runtime gate before adding the 1 MiB dispatch threshold.
- `SHA256SUMS`: hashes of public JSON/CSV/text artifacts in this directory.

Timing/profile JSONs are labeled public derivatives. They retain the original
report SHA-256 but replace the host label and reduce process observations to
counts and memory totals. The accepted comparison contains 200 process
observations, with no other compute process observed; this is not continuous
exclusive-device monitoring. Raw reports, result RDS files, failed attempts,
library checksums and raw CUDA traces are retained privately in
`validation/cell-sparse-pca-2026-10-19/benchmark-evidence-20261019.tar.gz`
(SHA-256 `a04631327bbb0943b5ef88feaefbd9f4f18ddf3dd8fcd27cf0677b85bb1ee8d2`).

The runtime passed 37 A100 assertions with native CUDA required and 36 local
native assertions. Its full Windows package check passed with 233 CPU-side
assertions, zero failures/warnings, and one expected native-hardware skip
inside R CMD check; the separate native runs cover that skipped test. Manuals
and vignettes passed. The check reported development/new-submission and one
5.83-second SingleCellExperiment example NOTE. Hosted platform checks remain
a separate integration gate and are recorded on the pull request.

The [performance guide](../../vignettes/performance-memory.Rmd) explains the
measured configurations, memory headroom, batch-size guidance, and next
profiling targets. The public 10x Genomics PBMC data is CC BY 4.0; its source
page and attribution are included there and in the input manifest.

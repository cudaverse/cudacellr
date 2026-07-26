# Changelog

## cudacellr 0.3.0

- Added
  [`cudacell_sce()`](https://cudaverse.github.io/cudacellr/reference/cudacell_sce.md)
  for a native, non-destructive `SingleCellExperiment` workflow.
  Normalized expression, PCA, HVG statistics, kNN relationships, and
  compute provenance are written to the corresponding assay, reduced
  dimension, `rowData`, `colPair`, and metadata locations.
- Existing assays, row and column metadata, reduced dimensions,
  alternative experiments, pairings, canonical size factors, labels, and
  user metadata are preserved. Output collisions fail before computation
  unless explicitly replaced with `overwrite = TRUE`.
- Delayed assays are never materialized silently. `realize = TRUE` is
  required to opt into an in-memory sparse realization.
- [`cuda_provenance()`](https://cudaverse.github.io/cudacellr/reference/cuda_provenance.md)
  now reads the ordered computation record directly from a
  `SingleCellExperiment` produced by
  [`cudacell_sce()`](https://cudaverse.github.io/cudacellr/reference/cudacell_sce.md).
- [`cuda_cell_pca()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_pca.md)
  and
  [`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md)
  now expose normalization `scale_factor`, `log1p`, and `min_mean`
  controls; the complete workflow also exposes PCA scaling through
  `scale.`.

## cudacellr 0.2.0

- Normalization, HVG selection, PCA, and kNN now compose into one
  ordered `cudaverse-stage/1` provenance record, including upstream
  sparse materialization.
- [`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md)
  reports its aggregate CPU/CUDA/hybrid compute device, and its print
  method no longer implies that the PCA device represents the complete
  workflow.
- Strict CUDA requests are validated before CPU preprocessing begins,
  avoiding expensive work before an unavailable-device error.
- Re-exported
  [`cuda_provenance()`](https://cudaverse.github.io/cudacellr/reference/cuda_provenance.md)
  as the common inspector for single-cell results.

## cudacellr 0.1.2

- Feature and cell identifiers now survive normalization, HVG selection,
  PCA, kNN, and the complete workflow. HVG results retain original
  feature indices, including when feature names are duplicated.
- Added a concise [`print()`](https://rdrr.io/r/base/print.html) method
  for workflow results and strict validation of the `log1p` flag.

## cudacellr 0.1.1

- [`cuda_cell_neighbors()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_neighbors.md)
  and
  [`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md)
  now expose `batch_size`, using the bounded-memory exact kNN
  implementation in `cudalearnr`.
- [`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md)
  now reuses its normalized matrix and variable-feature result for PCA
  instead of performing both preprocessing stages twice.

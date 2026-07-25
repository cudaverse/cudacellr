# Changelog

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

# Changelog

## cudacellr 0.1.1

- [`cuda_cell_neighbors()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_neighbors.md)
  and
  [`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md)
  now expose `batch_size`, using the bounded-memory exact kNN
  implementation in `cudalearnr`.
- [`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md)
  now reuses its normalized matrix and variable-feature result for PCA
  instead of performing both preprocessing stages twice.

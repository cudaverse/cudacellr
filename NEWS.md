# cudacellr 0.1.1

- `cuda_cell_neighbors()` and `cudacell_workflow()` now expose
  `batch_size`, using the bounded-memory exact kNN implementation in
  `cudalearnr`.
- `cudacell_workflow()` now reuses its normalized matrix and variable-feature
  result for PCA instead of performing both preprocessing stages twice.

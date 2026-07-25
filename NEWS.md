# cudacellr 0.1.2

- Feature and cell identifiers now survive normalization, HVG selection, PCA,
  kNN, and the complete workflow. HVG results retain original feature indices,
  including when feature names are duplicated.
- Added a concise `print()` method for workflow results and strict validation
  of the `log1p` flag.

# cudacellr 0.1.1

- `cuda_cell_neighbors()` and `cudacell_workflow()` now expose
  `batch_size`, using the bounded-memory exact kNN implementation in
  `cudalearnr`.
- `cudacell_workflow()` now reuses its normalized matrix and variable-feature
  result for PCA instead of performing both preprocessing stages twice.

# Package index

## Preprocess counts

- [`cuda_normalize_counts()`](https://cudaverse.github.io/cudacellr/reference/cuda_normalize_counts.md)
  : Normalize a single-cell count matrix
- [`cuda_hvg()`](https://cudaverse.github.io/cudacellr/reference/cuda_hvg.md)
  : Select highly variable features

## Reduce dimensions and find neighbours

- [`cuda_cell_pca()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_pca.md)
  : Run PCA on highly variable single-cell features
- [`cuda_cell_neighbors()`](https://cudaverse.github.io/cudacellr/reference/cuda_cell_neighbors.md)
  : Construct a nearest-neighbour representation

## End-to-end workflows

- [`cudacell_workflow()`](https://cudaverse.github.io/cudacellr/reference/cudacell_workflow.md)
  : Run the initial cudacellr workflow

## SingleCellExperiment integration

- [`cudacell_sce()`](https://cudaverse.github.io/cudacellr/reference/cudacell_sce.md)
  : Run a cudacell workflow on a SingleCellExperiment

## SeuratObject v5 integration

- [`cudacell_seurat()`](https://cudaverse.github.io/cudacellr/reference/cudacell_seurat.md)
  : Run a cudacell workflow on a Seurat object

## Provenance

- [`reexports`](https://cudaverse.github.io/cudacellr/reference/reexports.md)
  [`cuda_provenance`](https://cudaverse.github.io/cudacellr/reference/reexports.md)
  : Objects exported from other packages
- [`cuda_provenance(`*`<SingleCellExperiment>`*`)`](https://cudaverse.github.io/cudacellr/reference/cuda_provenance-methods.md)
  [`cuda_provenance(`*`<Seurat>`*`)`](https://cudaverse.github.io/cudacellr/reference/cuda_provenance-methods.md)
  : Inspect provenance stored on native single-cell containers

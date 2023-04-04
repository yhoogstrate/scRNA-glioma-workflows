#!/usr/bin/env R

check_config <- function(seurat_obj) {
  stopifnot(!is.null(seurat_obj@misc$dataset))
  stopifnot(!is.null(seurat_obj@misc$sample_name))
  stopifnot(!is.null(seurat_obj@misc$sample_full_name))
  stopifnot(!is.null(seurat_obj@misc$tumor_type))
  
  stopifnot(!is.null(seurat_obj@misc$max.percentage.mitochondrial))
  
  stopifnot(seurat_obj@misc$random_seed == 42)
  
  
  
  stopifnot(!is.null(seurat_obj$celltype_annotated))
  
  if(seurat_obj@misc$tumor_type == 'Astrocytoma, IDH-mut') {
    stopifnot(!is.null(seurat_obj$celltype_annotated_venteicher))
  }
}


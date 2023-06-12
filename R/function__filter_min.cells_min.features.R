#!/usr/bin/env R

filter_min_cells_min_features <- function(seurat_obj, min.cells = 0, min.features = 0) {
  
  if(min.features > 0) {
    nfeatures <- Matrix::colSums(x = Seurat::GetAssayData(object = seurat_obj, slot = "counts") > 0)
    seurat_obj <- seurat_obj[,which(x = nfeatures >= min.features)]
  }
  
  if(min.cells > 0) {
    num.cells <- Matrix::rowSums(x = Seurat::GetAssayData(object = seurat_obj, slot = "counts") > 0)
    seurat_obj <- seurat_obj[which(x = num.cells >= min.cells),]
  }
  
  return(seurat_obj)
}


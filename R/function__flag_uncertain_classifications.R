#!/usr/bin/env

flag_uncertain_classifications <- function(seurat_object, classifier='PSLDA') {
  n.excl <- 0
  
  slots <- data.frame(slot = names(seurat_object@misc$reclass.PSLDA.min.scores)) |> 
    dplyr::pull(slot)
  
  for(slot in slots) {
    key.slot <- paste0('reclass.',classifier,'.',slot)
    key.class <- paste0('reclass.',classifier,'.class')
    
    if(key.slot %in% names(seurat_object@meta.data)) {
      scores <- seurat_object@meta.data[[key.slot]]
      labels <- seurat_object@meta.data[[key.class]]
      
      min.score <- seurat_object@misc$reclass.PSLDA.min.scores[slot]
      
      exclude <- ifelse(is.na(labels), F,  (labels == slot) & (scores < min.score))
      
      seurat_object@meta.data[[key.class]] <- ifelse(exclude, NA, seurat_object@meta.data[[key.class]])
      
      n.excl <- n.excl + sum(exclude)
    }
  }
  
  print(paste0("excluded n=", n.excl, " ",classifier,  " prediction labels"))
  
  return (seurat_object)
}




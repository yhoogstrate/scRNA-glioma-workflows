#!/usr/bin/env R


geneset_component_association <- function(object, features, k.components) {
  
  df <- data.frame(PC = 1:k.components) |> 
    dplyr::mutate(p.val = pbapply::pblapply(PC, function(pc){
      df <- data.frame(PC = object@reductions$pca@feature.loadings[,pc]) |> 
        tibble::rownames_to_column('gene_symbol') |> 
        dplyr::mutate(col = gene_symbol %in% features)
      
      tt <- t.test(df[df$col == F,]$PC, df[df$col == T,]$PC , var.equal=T, mu=0)
      
      return(tt$p.value)
    } )) |> 
    dplyr::mutate(p.val = unlist(p.val)) |> 
    dplyr::mutate(PC.name = paste0("PC",PC))
  
  return(df)
}

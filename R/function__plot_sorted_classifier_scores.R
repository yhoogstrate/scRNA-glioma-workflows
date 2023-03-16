#!/usr/bin/env R


plot_sorted_classifier_scores <- function(seurat_object, col = "celltype.annotated", selector = "^reclass\\.PSLDA\\.") {
  
  slots <- data.frame(slot = names(seurat_object@meta.data)) |> 
    dplyr::filter(slot  == col | grepl(selector, slot))
  
  plt <- seurat_object@meta.data |>
    dplyr::select(slots$slot) |> 
    tidyr::pivot_longer(cols = -col) |> 
    dplyr::group_by(name) |> 
    dplyr::mutate(x = order(order(value))) |> 
    dplyr::ungroup()
  
  
  gg <- ggplot(plt, aes(x = x, y=value, col=celltype.annotated)) +
    geom_point(pch=19,cex=0.2) +
    facet_grid(rows = vars(name), scales="free",space="free") +
    labs(x=NULL, y=paste0(gsub("[^a-zA-Z0-9]+"," ",selector)," score / prob"))

  gg
  
  return(gg)

}



plot_sorted_classifier_scores(obj__diaz_2019__SF11136,'celltype.annotated' ,"^reclass\\.PSLDA\\.")

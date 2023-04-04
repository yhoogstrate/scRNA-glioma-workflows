#!/usr/bin/env R


plot_high_mitochondrial_cells <- function(seurat_obj) {
  
  plt <- data.frame(percentage_mitochondrial_reads = seurat_obj$percentage_mitochondrial_reads) |> 
    dplyr::mutate(filter = ifelse(percentage_mitochondrial_reads > seurat_obj@misc$max.percentage.mitochondrial, "excluded","included")) 
  
  
  ggobj <- ggplot(plt, aes(x = 1, y=percentage_mitochondrial_reads, col=filter)) +
    ggbeeswarm::geom_quasirandom(size=0.2, varwidth=T)  +
    geom_hline(yintercept = seurat_obj@misc$max.percentage.mitochondrial, lty=2, lwd=0.6) +
    ylim(0, 10)
    

  
  return (ggobj)
}


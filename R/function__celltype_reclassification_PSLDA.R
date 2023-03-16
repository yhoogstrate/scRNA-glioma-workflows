#!/usr/bin/env R


celltype_reclassification_pslda <- function(seurat_object, marker_features, n.components, slot_labels = "celltype.annotated", prefix="reclass.PSLDA.") {
  out <- data.frame()
  
  labels <- seurat_object[[slot_labels]] |> 
    tibble::rownames_to_column('umi')
  
  features <- seurat_object[marker_features,]@assays$RNA@data |> 
    as.matrix() |> 
    t() |> 
    as.data.frame(stringsAsFactors=F)

  stopifnot(nrow(labels) == nrow(features))
  stopifnot(labels$umi == features$umi)

  for(i in 1:10) {
    train.sel <- 1:nrow(data) %% 10 != (i - 1) 
    test.sel <- train.sel == F
    stopifnot(sum(train.sel) + sum(test.sel) == nrow(labels))
    
    
    train <- features[train.sel, ]
      #tibble::remove_rownames()
    
    train.labels <- labels[[slot_labels]][train.sel]
    stopifnot(nrow(train) == sum(train.sel))
    stopifnot(length(train.labels) == sum(train.sel))
    
    test <- features[test.sel,]
    stopifnot(nrow(test) == sum(test.sel))
    
    fit.plsda.caret <- caret::plsda(
      train,
      factor(train.labels),
      ncomp=n.components
      #,probMethod = "Bayes"
    )
    
    p.plsda.caret = caret:::predict.plsda(fit.plsda.caret, test, type = "prob") |> 
      as.data.frame() |>
      dplyr::rename_with( ~ gsub("\\.[0-9]+ comps$", "", .x, fixed = FALSE)) |> 
      tibble::rownames_to_column('umi') |> 
      (function(.) {
        assertthat::assert_that(nrow(.) == nrow(test))
        return(.)
      })()
  
    out <- rbind(out, p.plsda.caret)
  }
  
  
  stopifnot(nrow(data) == nrow(out))
  
  
  # reorder to original order and add prefix to colnames
  out <- data.frame(umi = rownames(data)) |> 
    dplyr::left_join(out, by=c('umi'='umi')) |> 
    tibble::column_to_rownames('umi') |> 
    dplyr::rename_with(~ paste0(prefix, .x))
  
  
  
  for(c in colnames(out)) {
    seurat_object[[c]] <- out[[c]]
  }
  
  
  return (seurat_object)
}



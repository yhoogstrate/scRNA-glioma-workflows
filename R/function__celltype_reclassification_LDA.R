#!/usr/bin/env R

celltype_reclassification_lda <- function(seurat_object, marker_features, slot_labels = "celltype.annotated", prefix="reclass.LDA.") {
  out <- data.frame()
  
  labels <- seurat_object[[slot_labels]] |> 
    tibble::rownames_to_column('umi')
  
  features <- seurat_object[marker_features,]@assays$RNA@data |> 
    as.matrix() |> 
    t() |> 
    as.data.frame(stringsAsFactors=F)

  stopifnot(nrow(labels) == nrow(features))
  stopifnot(labels$umi == features$umi)
  
  data <- features |> 
    tibble::rownames_to_column('umi') |> 
    dplyr::left_join(labels, by=c('umi'='umi'), suffix=c('','')) |> 
    tibble::column_to_rownames('umi')
  
  rm(labels, features)
  

  
  
  
  for(i in 1:10) {
    train.sel <- 1:nrow(data) %% 10 != (i - 1) 
    test.sel <- train.sel == F
    stopifnot(sum(train.sel) + sum(test.sel) == nrow(labels))
    
    
    train <- data[train.sel, ] |>
      dplyr::filter(get(slot_labels) != "NA") |>
      tibble::remove_rownames()
    stopifnot(nrow(train) == sum(train.sel))
    
    test <- data[test.sel,]
    stopifnot(nrow(test) == sum(test.sel))
    
    
    fit.lda <- MASS::lda(celltype.annotated ~ ., data=train)
    
    p.lda <- predict(
      fit.lda,
      test |>
        dplyr::select(-c(celltype.annotated))
    )
    
    add <- p.lda$posterior |> as.data.frame() |>  tibble::rownames_to_column('umi')
    stopifnot(nrow(add) == nrow(test))
    
    out <- rbind(out, add)
  }
  
  
  stopifnot(nrow(data) == nrow(out))
  
  
  # reorder to original order and add prefix to colnames
  out <- data.frame(umi = rownames(data)) |> 
    dplyr::left_join(out, by=c('umi'='umi')) |> 
    tibble::column_to_rownames('umi') |> 
    dplyr::rename_with( ~ paste0(prefix, .x))
  
  
  for(c in colnames(out)) {
    seurat_object[[c]] <- out[[c]]
  }
  
  
  return (seurat_object)
}




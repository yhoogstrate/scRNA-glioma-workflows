#!/usr/bin/env R

cache_seurat_object <- function(object) {
  if(is.na(object@misc$dataset)) {
    error("dataset slot not found in config")
  }

  if(is.na(object@misc$sample_name)) {
    error("sample_name slot not found in config")
  }

  if(is.na(object@misc$tumor_type)) {
    error("tumor_type slot not found in config")
  } else {
    tumor_type <- dplyr::recode(object@misc$tumor_type,
                                `Astrocytoma, IDH-mut` = 'A-IDH',
                                `Oligodendroglioma` = 'O-IDH',
                                `Glioblastoma` = 'GBM')
  }

  # ensure all code was committed
  if(length(system("git status --short", intern = TRUE)) == 0) {
    object@misc$git_branch <- system("git branch --show-current", intern = TRUE)
    object@misc$git_revision <- system("git rev-parse --short HEAD", intern = TRUE)
    
    fn <- paste0("cache/",
                 object@misc$dataset,
                 "__",
                 object@misc$sample_name,
                 "__[",
                 tumor_type,
                 "]__",
                 object@misc$git_branch,
                 ".Rds")
    
    saveRDS(object, file=fn)
    
    object@misc$git_branch <- NULL
    object@misc$git_revision <- NULL
    
  } else {
    warning("Not updating the cache object because not all sued code is commited to git")
  }
}

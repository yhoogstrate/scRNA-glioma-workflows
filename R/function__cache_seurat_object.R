#!/usr/bin/env R

cache_seurat_object <- function(object) {
  
  # ensure all code was committed
  if(length(system("git status --short", intern = TRUE)) == 0) {
    object@misc$git_branch <- system("git branch --show-current", intern = TRUE)
    object@misc$git_revision <- system("git rev-parse --short HEAD", intern = TRUE)

    fn <- paste0("cache/",
                 object@misc$sample_name,
                 "__",
                 object@misc$git_branch,
                 ".Rds")
    
    saveRDS(object, file=fn)
    
    object@misc$git_branch <- NULL
    object@misc$git_revision <- NULL
    
  } else {
    warning("Not updating the cache object because not all sued code is commited to git")
  }
}

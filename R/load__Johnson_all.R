
# install libs ----

## Seurat ----

## SeuratData ----

#devtools::install_github('satijalab/seurat-data')

## SeuratDisk ----

#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#remotes::install_github("mojaveazure/seurat-disk")


# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html



if(!file.exists("cache/analysis_scRNAseq_tumor_counts.h5")) {
  SeuratDisk::Convert("data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_counts.h5ad", dest = "cache/analysis_scRNAseq_tumor_counts.h5", overwrite = TRUE)
}
seurat_obj <- SeuratDisk::LoadH5Seurat("cache/analysis_scRNAseq_tumor_counts.h5")

rename <- read.table('cache/syn25956426_Johnson_cell_identifiers.txt', header=T) |> 
  dplyr::rename(new.name = umi) |> 
  dplyr::mutate(old.name = paste0("Cell",1:dplyr::n())) 
seurat_obj <- Seurat::RenameCells(seurat_obj, new.names = rename$new.name)

metadata <- rename |> 
  dplyr::left_join(read.csv("data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_metadata.tsv", sep="\t"), by=c('new.name'='cell_barcode')) |> 
  dplyr::left_join(
    read.csv('data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_syn25880693_clinical_metadata.csv'),
    by=c('case_barcode'='case_barcode')
  )



stopifnot(colnames(seurat_obj) == metadata$new.name)

for(slot in c("cell_state","case_barcode","case_sex","tumor_location","laterality","driver_mutations","time_point","idh_codel_subtype","who_grade","histological_classification")) {
  seurat_obj[[slot]] <- metadata[[slot]]
}



# export individual samples ----

# split per individual and export
for(slot in read.csv('data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_syn25880693_clinical_metadata.csv')$case_barcode) {
  #slot = "SM001"
  
  export <- subset(seurat_obj, case_barcode == slot)
  export.slim <- Seurat::DietSeurat(export, 
                                    counts = TRUE,
                                    data = TRUE,
                                    scale.data = FALSE,
                                    features = NULL,
                                    assays=c('RNA'),
                                    dimreducs = NULL,
                                    graphs = NULL,
                                    misc = F
                                    )
  
  export.slim@misc <- export.slim@meta.data[c('case_barcode', 'case_sex', 'tumor_location', 'laterality', 'driver_mutations', 'time_point', 'idh_codel_subtype', 'who_grade')] |> 
    tibble::remove_rownames() |> 
    dplyr::distinct()
  
  export.slim@meta.data <- export.slim@meta.data[c('cell_state')]
  
  DropletUtils::write10xCounts(path=paste0("cache/syn25956426_Johnson__split__",slot), x=export.slim@assays$RNA@counts)
  
  write.csv(data.frame(
    cell_barcode = colnames(export.slim),
    cell_state = export.slim@meta.data$cell_state
  ), file=paste0("cache/syn25956426_Johnson__split__",slot,"/cell_states.csv"))
  
  rm(export, export.slim)
}





# confirm validaty labels ----
# 
# # norm
seurat_obj <- Seurat::NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = F)
seurat_obj <- Seurat::FindVariableFeatures(object = seurat_obj, selection.method = "vst", verbose = F)
seurat_obj <- Seurat::ScaleData(object = seurat_obj, verbose = T)
seurat_obj <- Seurat::RunPCA(reduction.key = "PC_", object = seurat_obj, features = Seurat::VariableFeatures(object = seurat_obj),verbose = F)
seurat_obj <- Seurat::FindNeighbors(object = seurat_obj, dims = 1:35, verbose = F)
seurat_obj <- Seurat::FindClusters(object = seurat_obj, resolution = 1.2, algorithm = 1, verbose = F)
seurat_obj <- Seurat::RunUMAP( object = seurat_obj, dims = 1:35, verbose = T )


# Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_state")
# Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .6, group.by = "case_barcode")
# Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .6, group.by = "idh_codel_subtype")


Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_state")

# Stem-like = tumor

Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_state")

# Diff.-like = tumor
# Prolif. stem-like = tumor

# Oligodendrocyte


# Myeloid
# - Dendritic cell
# Granulocyte
Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_state")
Seurat::FeaturePlot(seurat_obj, features=c('CD163' ,'P2RY12' ,'CD14' ,'ITGB2' ,'C1QC' ,'SLC11A1'))

# T cell - klopt

Seurat::FeaturePlot(seurat_obj, features=c('CD2','CD3D','TRBC2','TRAC'))

# B cell
Seurat::FeaturePlot(seurat_obj, features = c("IGLC3","CD19","CD79B"),order=T)



## endothelial & pericytes
Seurat::FeaturePlot(seurat_obj,
                    reduction = "umap",
                    features = celltype_markers |>
                      dplyr::filter(celltype == "endothelial") |> 
                      dplyr::pull(markers) |> 
                      head(n = 4),
                    label=F)


Seurat::FeaturePlot(seurat_obj,
                    reduction = "umap",
                    features = celltype_markers |>
                      dplyr::filter(celltype == "pericyte") |> 
                      dplyr::pull(markers) |> 
                      head(n = 4),
                    label=F)


Idents(seurat_obj)<- "cell_state"
m.fibroblast2 <- Seurat::FindMarkers(seurat_obj, ident.1 = 'Fibroblast', ident.2='Pericyte')#, slot="cell_state"
Idents(seurat_obj)<- "seurat_clusters"


Seurat::FeaturePlot(seurat_obj,
                    reduction = "umap",
                    features = c('CYP1B1','OGN','CLU','CXCL12'),
                    label=F)






# Fibroblast
# Pericyte
# Endothelial
Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_state")



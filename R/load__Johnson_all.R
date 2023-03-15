
# install libs ----

## Seurat ----

## SeuratData ----

devtools::install_github('satijalab/seurat-data')

## SeuratDisk ----

#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
remotes::install_github("mojaveazure/seurat-disk")


# load libs

library(Seurat)
library(SeuratData)
library(SeuratDisk)


# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html


SeuratDisk::Convert("data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_counts.h5ad", dest = "cache/analysis_scRNAseq_tumor_counts.h5", overwrite = TRUE)
d = SeuratDisk::LoadH5Seurat("cache/analysis_scRNAseq_tumor_counts.h5")
colnames(d)[c(1:5,55248)]
d <- Seurat::NormalizeData( object = d,  normalization.method = "LogNormalize", scale.factor = 1e4, verbose = F)





m = read.csv("data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_metadata.tsv",sep="\t")
c = read.csv("data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_gene_expression_20211001_full.tsv",sep="\t")

head(c)


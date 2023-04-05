
#expression_matrix_full <- read.delim("data/output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#")
#Load gene annotation data
if(!file.exists('cache/gencode.31.Rds')) {
  gencode.31 <- read.delim("data/ref/star-hg19/gencode.v31lift37.annotation.gtf", comment.char="#",stringsAsFactors = F,header=F) |> 
    dplyr::filter(V3 == "gene") |> 
    dplyr::mutate(ENSG = gsub("^.+(ENSG[^;]+);.+$","\\1",V9)) |> 
    dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1",V9)) |> 
    dplyr::mutate(gene_type = gsub("^.+gene_type ([^;]+);.+$","\\1",V9)) |> 
    dplyr::mutate(V9 = NULL)
  
  saveRDS(gencode.31, 'cache/gencode.31.Rds')
} else {
  gencode.31 <- readRDS('cache/gencode.31.Rds')
}






sid <- 'L_CPT0205890014'
origin_file <-  "data/CPTAC-3/L_CPT0205890014/2f54e97e-0961-410a-b429-f42d212c4755/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}
rm(origin_file, sid)

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version


mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- Seurat::CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)
rm(mat, lfile, gene.order, barcode, tt)


saveRDS(object_1, file="cache/CPTAC-3__L_CPT0205890014_C3N-02783_export.Rds")





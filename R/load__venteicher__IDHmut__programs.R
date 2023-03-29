#!/usr/bin/env R

# From: Decoupling genetics, lineages, and microenvironment in IDH-mutant gliomas by single-cell RNA-seq
# https://pubmed.ncbi.nlm.nih.gov/28360267/

celltype_markers_Venteicher_IDHmut <- read.table("assets/celltype_markers_Venteicher_IDHmut.txt",header=T) |> 
  dplyr::rename(gene_symbol = Gene) |> 
  assertr::verify(is.logical(Oligo.program)) |> 
  assertr::verify(is.logical(Astro.program)) |> 
  assertr::verify(is.logical(Stemness.program))


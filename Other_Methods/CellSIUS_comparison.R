library(Seurat)
library(Ckmeans.1d.dp)
library(data.table)
library(igraph)
library(ggplot2)
library(Rcpp)



setwd("/home/ies/gabriele.lubatti/entropy_mixing")
load("elmir_seurat.Rda")
umap_elmir=readRDS("annot_umap.rds")
norm_elmir=as.matrix(GetAssayData(elmir_seurat, slot = "data",assay="RNA"))


source("/home/ies/gabriele.lubatti/entropy_mixing/cellsius/ane/CellSIUS.R")

CellSIUS_elmir_new<-CellSIUS(mat.norm = norm_elmir,group_id = umap_elmir$cluster_id,min_n_cells=5, min_fc = 1,mcl_path = "/home/ies/gabriele.lubatti/local/bin/mcl")

setwd("/home/ies/gabriele.lubatti/entropy_mixing/cellsius/ane")
save(CellSIUS_elmir_new,file="CellSIUS_elmir_new.Rda")
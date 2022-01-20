

Cluster_analysis_integrate_rare=function(raw_counts,batch_name_1,resolution,neighbors,max_dimension,feature_genes=NULL){
  batch_1 <- CreateSeuratObject(counts = raw_counts, project = batch_name_1)
  batch_1 <- subset(batch_1, subset = nFeature_RNA > 0)
  batch_1 <- NormalizeData(batch_1, verbose = FALSE)
  batch_1 <- FindVariableFeatures(batch_1, selection.method = "vst", nfeatures = 2000)

  batch_combined=batch_1



  batch_combined <- ScaleData(batch_combined, verbose = FALSE)

  batch_combined <- RunPCA(batch_combined, npcs = max_dimension, verbose = FALSE,features=feature_genes)
  batch_combined=RunUMAP(batch_combined, reduction = "pca", dims = 1:20,set.seed=42)
  batch_combined <- FindNeighbors(batch_combined, reduction = "pca", dims = 1:max_dimension,k.param=neighbors)
  batch_combined <- FindClusters(batch_combined, resolution =resolution)
  return(batch_combined)}



find_resolution=function(seurat_object,resolution_vector){
  for (i in resolution_vector){
    seurat_object=FindClusters(seurat_object, resolution =i)
  }
  clustree(seurat_object)
}






Function_fisher=function(a,b,c){
  matrice=matrix(c(sum(c[c%in%a]%in%c[c%in%b]),sum(c[c%in%a]%in%c[!c%in%b]),sum(c[!c%in%a]%in%c[c%in%b]),sum(c[!c%in%a]%in%c[!c%in%b])),nrow=2,ncol=2,byrow = T)
  fisher.test(matrice)
}






De_seurat_cluster=function(seurat_object,cluster,names_cell,max_p_value){
  level=levels(as.factor(cluster))
  final_markers=vector("list",length(level))
  for (i in 1:length(level)){
    markers=FindMarkers(seurat_object,ident.1=names_cell[cluster==level[i]],ident.2=names_cell[cluster!=level[i]],only.pos = T)
    markers_new=markers[markers$p_val_adj<=max_p_value,]
    markers_new=markers_new[order(markers_new$p_val_adj),]
    markers_final=row.names(markers_new)
    final_markers[[i]]=markers_final

  }

  return(final_markers)
}









merge_cluster=function(old_cluster,new_cluster,max_number=NULL){

  if(is.null(max_number)){
    cluster_final=old_cluster
    cluster_final[match(names(new_cluster),names(old_cluster))]=new_cluster
    return(cluster_final)
  }
  else{
    index_small=which(table(new_cluster)<max_number)
    cluster_small=names(table(new_cluster))[index_small]
    print(cluster_small)
    cluster_final=old_cluster
    for (i in 1:length(cluster_small)){
      cluster_final[names(old_cluster)%in%names(new_cluster)[as.vector(new_cluster)==cluster_small[i]]]=paste(cluster_small[i],"step_2",sep="-")
    }
    return(cluster_final)
  }
}



test_hvg=function (raw_counts, final_cluster_version, relevant_genes,
                   background, number_hvg,min_p_value)
{
  levels_cluster = levels(as.factor(final_cluster_version))
  final_p_value = rep(list(0), length(levels))
  for (i in 1:length(levels_cluster)) {
    raw_counts_small = raw_counts[, final_cluster_version ==
                                    levels_cluster[i]]
    batch_1 <- CreateSeuratObject(counts = raw_counts_small,
                                  project = "sub_cluster")
    batch_1 <- subset(batch_1, subset = nFeature_RNA > 0)
    batch_1 <- NormalizeData(batch_1, verbose = FALSE)
    batch_1 <- FindVariableFeatures(batch_1, selection.method = "vst",
                                    nfeatures = number_hvg)
    hvg_genes = batch_1@assays$RNA@var.features
    fisher_output = (Function_fisher(relevant_genes, hvg_genes,
                                     background))
    final_p_value[[i]] = list(fisher_output[[1]], fisher_output[[3]],
                              relevant_genes[relevant_genes %in% hvg_genes])
  }
  names(final_p_value) = levels_cluster
  select_sub=lapply(final_p_value, function(x){if(x[[1]]<min_p_value & x[[2]]>1){
    return(TRUE)}
    else{return(FALSE)}
  })

  final_sub=names(select_sub)[select_sub==TRUE]
  return(list(final_p_value,final_sub))
}














Cluster_analysis_sub=function (raw_counts, resolution, neighbors, max_dimension,
                                    name_cluster)
{
  batch_1 <- CreateSeuratObject(counts = raw_counts, project = "sub_cluster")
  batch_1 <- subset(batch_1, subset = nFeature_RNA > 0)
  batch_1 <- NormalizeData(batch_1, verbose = FALSE)
  batch_1 <- FindVariableFeatures(batch_1, selection.method = "vst",
                                  nfeatures = 2000)
  batch_combined <- batch_1
  batch_combined <- ScaleData(batch_combined, verbose = FALSE)
  hvg_genes = batch_1@assays$RNA@var.features
  #genes_cluster = intersect(relevant_genes, hvg_genes)
  #print(length(genes_cluster))
  genes_cluster=hvg_genes
  batch_combined <- RunPCA(batch_combined, npcs = max_dimension,
                           verbose = FALSE, features = genes_cluster)
  #batch_combined = RunUMAP(batch_combined, reduction = "pca",
  #dims = 1:20, set.seed = 42)
  batch_combined <- FindNeighbors(batch_combined, reduction = "pca",
                                  dims = 1:max_dimension, k.param = neighbors)
  batch_combined <- FindClusters(batch_combined, resolution = resolution)
  new_cluster = as.vector(batch_combined$seurat_clusters)
  for (i in 1:length(levels(batch_combined$seurat_clusters))) {
    new_cluster[as.vector(batch_combined$seurat_clusters) ==
                  levels(batch_combined$seurat_clusters)[i]] = paste(name_cluster,
                                                                     new_cluster[as.vector(batch_combined$seurat_clusters) ==
                                                                                   levels(batch_combined$seurat_clusters)[i]], sep = "_")
  }
  names(new_cluster) = names(batch_combined$seurat_clusters)
  batch_combined$seurat_clusters = new_cluster
  return(batch_combined)
}





markers_first_method=function(seurat_object,clusters,cell_names,number_top){
  level=levels(as.factor(clusters))
  marker_gabriele=as.list(rep(0,length(level)))
  marker_all=De_seurat_cluster(seurat_object,clusters,cell_names,0.05)

  for (i in 1:length(level)){

    marker_gabriele[[i]]=marker_all[[i]]




    print(i)

  }
  marker_all=unlist(marker_gabriele)
  marker_all=marker_all[isUnique(marker_all)]

  for (i in 1:length(level)){
    marker_gabriele[[i]]=marker_gabriele[[i]][marker_gabriele[[i]]%in%marker_all]
  }


  marker_to_see=as.list(rep(0,length(level)))
  marker_complete=as.list(rep(0,length(level)))
  for (i in 1:length(level)){
    if (length(marker_gabriele[[i]])>=number_top){
      marker_to_see[[i]]=marker_gabriele[[i]][1:number_top]
      names(marker_to_see[[i]])=rep(level[i],length(marker_to_see[[i]]))
      marker_complete[[i]]=marker_gabriele[[i]]
      names(marker_complete[[i]])=rep(level[i],length(marker_complete[[i]]))
    }
    if ((length(marker_gabriele[[i]])>0)&(length(marker_gabriele[[i]])<number_top)){
      marker_to_see[[i]]=marker_gabriele[[i]][1:length(marker_gabriele[[i]])]
      names(marker_to_see[[i]])=rep(level[i],length(marker_to_see[[i]]))
      marker_complete[[i]]=marker_gabriele[[i]]
      names(marker_complete[[i]])=rep(level[i],length(marker_complete[[i]]))
    }
    if (length(marker_gabriele[[i]])==0){
      marker_to_see[[i]]=NULL
      marker_complete[[i]]=NULL
    }
  }



  marker_plot=unlist(marker_to_see)
  marker_plot=marker_plot[marker_plot!=0]
  marker_complete=unlist(marker_complete)
  marker_complete=marker_complete[marker_complete!=0]
  marker_plot_plot_sing=as.list(rep(0,length(level)))


  for (i in 1:length(level)){
    if(sum(names(marker_plot)==level[i])>0){
      marker_plot_plot_sing[[i]]=paste(marker_plot[names(marker_plot)==level[i]],level[i],sep="_")
    }
    else{}
  }
  marker_plot_plot=unlist(marker_plot_plot_sing)
  marker_plot_plot=marker_plot_plot[marker_plot_plot!=0]

  return(list(marker_plot,marker_plot_plot,marker_complete))}





white_black_marker=function(cluster_race_entropy,single_cluster,norm_race,marker_list,threshold=0){

  white_black_7=apply(norm_race[marker_list[names(marker_list)==single_cluster],],1,function(x){
    mean_one=median(x[cluster_race_entropy==single_cluster])
    mean_different=rep(0,length(levels(factor(cluster_race_entropy[cluster_race_entropy!=single_cluster]))))
    for ( i in 1:length(levels(factor(cluster_race_entropy[cluster_race_entropy!=single_cluster])))){
      mean_different[i]=median(x[cluster_race_entropy==levels(factor(cluster_race_entropy[cluster_race_entropy!=single_cluster]))[i]])}

    if (mean_one>threshold&sum(mean_different)==0){
      return(TRUE)
    }
    else{return(FALSE)}
  })
  return(white_black_7)
}


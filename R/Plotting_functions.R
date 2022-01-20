













Plot_umap=function(cordinate_umap,cluster){


  umap_plot=ggplot(cordinate_umap, aes(cordinate_umap[,1],cordinate_umap[,2] )) +
    geom_point(aes(colour =as.factor(cluster)))+xlab("UMAP 1")+ylab("UMAP 2")+labs(col="Cluster")



  return(list(umap_plot))
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


Plot_genes_sum=function(cordinate_umap,norm_counts,genes_relevant,name_title){

  norm_counts_small=apply(norm_counts[genes_relevant,],1,function(x){
    y=x/sum(x)
    return(y)
  })


  gene_sum=apply(norm_counts_small,1,sum)



  index_sort=order(gene_sum)
  row.names(cordinate_umap)=colnames(norm_counts)
  cordinate_umap=cordinate_umap[colnames(norm_counts)[index_sort],]
  gene_sum=sort(gene_sum)

  umap_plot=ggplot(cordinate_umap, aes(cordinate_umap[,1],cordinate_umap[,2] )) +
    geom_point(aes(colour =(gene_sum)))+xlab("UMAP 1")+ylab("UMAP 2")+labs(col="Log norm counts")+
    scale_colour_gradient(low = "white", high = "blue")+
    ggtitle(name_title)
  return((umap_plot))
}


Plot_gene=function(norm_counts,cordinate_umap,gene_id,title_name){


  gene_expr=as.vector(norm_counts[gene_id,])
  index_sort=order(gene_expr)
  row.names(cordinate_umap)=colnames(norm_counts)
  cordinate_umap=cordinate_umap[colnames(norm_counts)[index_sort],]
  gene_expr=sort(gene_expr)

  umap_plot=ggplot(cordinate_umap, aes(cordinate_umap[,1],cordinate_umap[,2] )) +
    geom_point(aes(colour =(gene_expr)))+xlab("UMAP 1")+ylab("UMAP 2")+labs(col="Log norm counts")+
    scale_colour_gradient(low = "white", high = "blue")+
    ggtitle(title_name)
  return((umap_plot))
}







heatmap_marker_all_new=function(marker_plot,marker_plot_plot,my.clusters,Condition,norm_counts,number){

  cluster_unique=unique(my.clusters)
  cluster_unique=sort(cluster_unique,decreasing = F)

  index=as.list(rep(0,length(unique(my.clusters))))
  for (i in 1:length(cluster_unique)){
    index[[i]]=which(my.clusters==cluster_unique[i])
  }

  condition_unique=as.list(rep(0,length(unique(my.clusters))))
  for (i in 1:length(cluster_unique)){
    condition_unique[[i]]=Condition[index[[i]]]
  }

  cell_all_day=colnames(norm_counts)
  cell_unique=as.list(rep(0,length(unique(my.clusters))))
  for (i in 1:length(cluster_unique)){
    cell_unique[[i]]=cell_all_day[index[[i]]]
  }

  cluster_unique_2=as.list(rep(0,length(unique(my.clusters))))
  for (i in 1:length(cluster_unique)){
    cluster_unique_2[[i]]=my.clusters[index[[i]]]
  }

  cluster_finale=unlist(cluster_unique_2)

  good_finale=unlist(cell_unique)
  medium_finale=unlist(condition_unique)

  norm_counts_plot=norm_counts[marker_plot,good_finale]

  row.names(norm_counts_plot)<-marker_plot_plot

  color_cluster=rep(0,length(unique(my.clusters)))

  for (i in 1:length(color_cluster) ){
    color_cluster[i]=gg_color_hue(length(color_cluster))[i]
  }
  names(color_cluster)=as.character(cluster_unique)


  color_condition=rep(0,length(unique(Condition)))

  for (i in 1:length(color_condition) ){
    color_condition[i]=gg_color_hue(length(color_condition))[i]
  }
  names(color_condition)=as.character(unique(Condition))

  Condition_factor=factor(Condition)

  names(color_condition)=as.character(levels(Condition_factor))



  data_heatmap=data.frame(Cluster = cluster_finale
                          ,Condition = medium_finale)


  haCol1 = HeatmapAnnotation(df= data_heatmap,col=list(Cluster=color_cluster,Condition=color_condition)

                             , show_legend = T
  )

  ht21 = Heatmap(as.matrix(norm_counts_plot)# heat.vals_Mutant_Paired #
                 ,cluster_rows = FALSE
                 , col= colorRamp2(c(0, round(max(norm_counts_plot))), c("white", "red"))
                 , name = "log norm counts"
                 ,cluster_columns = F
                 , top_annotation = haCol1
                 , row_names_gp = gpar(fontsize = number
                 )
                 , show_column_names = F
                 , show_row_names = T
  )

  draw(ht21

  )



}









interactive_plot=function(coordinate_umap,Color,text,min_x=NULL,max_x=NULL,min_y=NULL,max_y=NULL){
  colnames(coordinate_umap)=c("UMAP_1","UMAP_2")
  fig <- plot_ly(data = coordinate_umap,
                 x = ~UMAP_1, y = ~UMAP_2,
                 color = ~Color,
                 type = "scatter",
                 mode = "markers",
                 marker = list(size = 5, width=2),
                 text=~text,
                 hoverinfo="text"
  )

  if(!is.null(min_x)){

    fig <- fig %>% layout(
      xaxis = list(range = c(min_x, max_x)),
      yaxis = list(range = c(min_y, max_y))
    )


  }

  return(fig)
}












selection_localized_genes=function (norm_counts, relevant_genes, min_number_cells = 4,
                                           max_number_genes = 10)
{
  norm_counts_entropy = norm_counts[(relevant_genes),
                                                ]
  logic_entropy = apply(norm_counts_entropy, 1, function(x) {
    y = as.list(x)
    return(y)
  })
  logic_entropy_greater = lapply(logic_entropy, function(x) {
    y = x[!all(x == 0)]
    if (length(y) > 0) {
      return(y)
    }
    else {
      return(NULL)
    }
  })
  index = rep(0, length(logic_entropy_greater))
  for (i in 1:length(logic_entropy_greater)) {
    index[i] = sum(logic_entropy_greater[[i]] > 0)
  }
  top_localized = (relevant_genes)
  top_localized = top_localized[index > min_number_cells]
  top_genes = apply(norm_counts[top_localized, ], 2,
                    function(x) {
                      y = sort(x, decreasing = T)
                      z = y[1:10]
                      k = top_localized[order(x, decreasing = T)][1:max_number_genes]
                      return(list(z, k))
                    })
  final_top_genes = top_genes
  somma = rep(0, length(final_top_genes))
  for (i in 1:length(final_top_genes)) {
    somma[i] = sum(as.numeric(final_top_genes[[i]][[1]]) >
                     1)
  }
  text = rep("NO", length(colnames(norm_counts)))
  for (i in 1:length(final_top_genes)) {
    text[i] = paste(final_top_genes[[i]][[2]], collapse = "-")
  }
  text[somma == 0] = "No specific genes"
  return(text)
}


Plot_balons=function (norm_counts_ane_24, final_cluster_version_sub, marker_plot_day_24_first_plot,
                          max_number, max_size = 5, text_size = 7)
{
  livelli = levels(factor(final_cluster_version_sub))
  all_markers_plot = rep(list(0), length(livelli))
  all_markers = rep(list(0), length(livelli))
  for (i in 1:length(livelli)) {
    genes_4_0 = marker_plot_day_24_first_plot[(names(marker_plot_day_24_first_plot) ==
                                                 livelli[i])]
    if (length(genes_4_0) > 0) {
      genes_4_0 = genes_4_0[1:max_number]
      all_markers[[i]] = genes_4_0
      genes_4_0_plot = paste(livelli[i], genes_4_0, sep = "_")
      all_markers_plot[[i]] = genes_4_0_plot
    }
  }
  all_markers = lapply(all_markers, function(x) {
    x[x != 0]
  })
  all_markers_plot = lapply(all_markers_plot, function(x) {
    x[x != 0]
  })
  all_markers_values = rep(list(0), length(livelli))
  for (i in 1:length(livelli)) {
    gene_values_0 = apply(norm_counts_ane_24[unlist(all_markers,
                                                    use.names = FALSE), final_cluster_version_sub ==
                                               livelli[i]], 1, mean)
    all_markers_values[[i]] = gene_values_0
  }
  values = rep(list(0), length(livelli))
  for (i in 1:length(livelli)) {
    valore = apply(norm_counts_ane_24[unlist(all_markers,
                                             use.names = FALSE), final_cluster_version_sub ==
                                        livelli[i]], 1, function(x) {
                                          sum(x > 1)/length(x)
                                        })
    values[[i]] = valore
  }
  cluster_label_all = rep(list(0), length(livelli))
  for (i in 1:length(livelli)) {
    cluster_label_all[[i]] = rep(livelli[i], length(unlist(all_markers)))
  }
  all_markers_plot = unlist(all_markers_plot)
  cluster_label_all = unlist(cluster_label_all)
  values = unlist(values)
  all_markers_values = unlist(all_markers_values)
  balloon_melted = data.frame(all_markers_plot, cluster_label_all,
                              values)
  colnames(balloon_melted) = c("genes_plot", "cluster_label_all",
                               "values")
  p <- ggplot(balloon_melted, aes(x = cluster_label_all, y = as.character(genes_plot)))
  p + geom_point(aes(size = values, colour = as.numeric(all_markers_values))) +
    theme(panel.background = element_blank(), panel.border = element_rect(colour = "blue",
                                                                          fill = NA, size = 3)) + scale_size_area(max_size = max_size) + scale_colour_gradient(low = "white",
                                                                                                                                                               high = "midnight blue", na.value = NA) + labs(colour = "Mean log2 norm expression",
                                                                                                                                                                                                             size = "Value (fraction cells above 1 log2 norm counts)") +
    xlab("Cluster") + ylab("Markers")+theme(text = element_text(size=text_size),axis.text.x = element_text(angle=45,vjust=1,hjust=1),axis.title.x=element_blank())
}

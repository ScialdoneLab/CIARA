


Get_background_full=function(norm_matrix,threshold=1,n_cells_low=3,n_cells_high=20){

  genes_filter=apply(norm_matrix,1,function(x){
    x=x[x>threshold]
    if (length(x)>=n_cells_low){return(TRUE)}
    else{return(FALSE)}
  })

  genes_filter_2=apply(norm_matrix,1,function(x){
    x=x[x>threshold]
    if (length(x)>n_cells_high){return(FALSE)}
    else{return(TRUE)}
  })
  genes_important=row.names(norm_matrix)[genes_filter]
  genes_important_2=row.names(norm_matrix)[genes_filter_2]
  genes_important_final=intersect(genes_important,genes_important_2)
  return(genes_important_final)
}








CIARA_gene=function(norm_matrix,knn_matrix,gene_expression,p_value=0.001,odds_ratio=2,local_region=1,approximation=FALSE){




  median_genes=median(gene_expression)
  binary_expression=rep(0,length(gene_expression))
  binary_expression[gene_expression>median_genes]=1
  binary_expression[gene_expression<=median_genes]=0


  if (approximation==FALSE){
    sub_feature=colnames(norm_matrix)} else {
      sub_feature=colnames(norm_matrix)[which(binary_expression==1)]
    }



  p_value_gene=rep(1,length(sub_feature))
  knn_matrix_big=knn_matrix[sub_feature,]



  perform_sub_entropy=function(j){

    nn_cell=colnames(norm_matrix)[which(knn_matrix_big[sub_feature[j],]>0)]
    nn_gene_expression=binary_expression[colnames(norm_matrix)%in%nn_cell]
    if(sum(nn_gene_expression)!=0){
      input_fisher=matrix(c(sum(nn_gene_expression==1),sum(binary_expression==1)-sum(nn_gene_expression==1),sum(nn_gene_expression==0),sum(binary_expression==0)-sum(nn_gene_expression==0)),nrow=2,ncol=2,byrow = T)
      test_output=fisher.test(input_fisher,alternative = "g")
      p_value_gene[j]=test_output$p.value



    }

    return((p_value_gene[j]))
  }

  p_value_gene=unlist(lapply(seq(1:length(sub_feature)),perform_sub_entropy))

  if(sum(p_value_gene<p_value)>=local_region){
    p_value_final=min(p_value_gene)}
  else{
    p_value_final=1
  }

  return((p_value_final))
}


CIARA=function(norm_matrix,knn_matrix,background,cores_number=1,p_value=0.001,odds_ratio=2,local_region=1,approximation=FALSE){
  run_loop_genes=function(i){
    gene_expression=as.vector(norm_matrix[background[i],])
    print(paste("Gene:",background[i],"number:",i,sep=" "))
    return(CIARA_gene(norm_matrix,knn_matrix,gene_expression,p_value,odds_ratio,local_region,approximation))

  }



  result_final=do.call(rbind, mclapply(seq(1:length(background)),run_loop_genes , mc.cores = cores_number))


  row.names(result_final)=background

  return(result_final)
}
























#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
gc()
rm(list=ls())
library(optparse)
option_list = list(
  make_option(c("-n", "--norm_matrix"), type="character", default=NULL,
              help="normalized count matrix", metavar="character"),
  make_option(c("-k", "--knn_matrix"), type="character", default=NULL,
              help="k nearest neighbors matrix", metavar="character"),
  make_option(c("-r", "--local_region"), type="integer", default=1,
              help="minimum number of local region with an enrichment for a given gene [default= %default]", metavar="integer"),
  make_option(c("-f", "--odds_ratio"), type="integer", default=2,
              help="Minimum odds ratio from fisher test for having enrichemnt [default= %default]", metavar="integer"),
  make_option(c("-l", "--n_cells_low"), type="integer", default=3,
              help="minimum number cells where a gene is expressed at a level above threshold  [default= %default]", metavar="integer"),
  make_option(c("-h", "--n_cells_high"), type="integer", default=20,
              help="maximum number of cells where a gene is expressed at a level above threshold  [default= %default]", metavar="integer"),
  make_option(c("-t", "--threshold"), type="integer", default=1,
              help="threshold in expression for a given gene [default= %default]", metavar="integer"),
  make_option(c("-o", "--output"), type="character", default="out.Rda",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-p", "--min_p_value"), type="double", default=0.001,
              help="Threshold for p value for fisher test[default= %default]", metavar="double"),
  make_option(c("-c", "--cores_number"), type="integer", default=10,
              help="Number of cores to use[default= %default]", metavar="integer"),
  make_option(c("-a", "--approximation"), type="logical", default=FALSE,
              help="For a given gene, the fisher test is run in the local regions of only the cells where the binarized gene expression is 1[default= %default]", metavar="logical")

)


opt_parser = OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$'norm_matrix')|is.null(opt$'knn_matrix') ){
  print_help(opt_parser)
  stop("Two arguments must be supplied (norm counts and k nearest neighbors ).n", call.=FALSE)
}



library(DescTools)
library(rlist)
library(parallel)












Get_background_full=function(norm_matrix,threshold,n_cells_low,n_cells_high){

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



CIARA_gene=function(norm_matrix,knn_matrix,gene_expression,p_value,odds_ratio,local_region,approximation){




  median_genes=median(gene_expression)
  binary_expression=rep(0,length(gene_expression))
  binary_expression[gene_expression>median_genes]=1
  binary_expression[gene_expression<=median_genes]=0


  if (approximation==FALSE){
    sub_feature=colnames(norm_matrix)
    } else {
      sub_feature=colnames(norm_matrix)[which(binary_expression==1)]
    }



  #p_value_gene=rep(1,length(sub_feature))
  knn_matrix_big=knn_matrix[sub_feature,]



  perform_sub_entropy=function(j){

    nn_cell=colnames(norm_matrix)[which(knn_matrix_big[sub_feature[j],]>0)]
    nn_gene_expression=binary_expression[colnames(norm_matrix)%in%nn_cell]
    if(sum(nn_gene_expression)!=0){
      input_fisher=matrix(c(sum(nn_gene_expression==1),sum(binary_expression==1)-sum(nn_gene_expression==1),sum(nn_gene_expression==0),sum(binary_expression==0)-sum(nn_gene_expression==0)),nrow=2,ncol=2,byrow = T)
      test_output=fisher.test(input_fisher,alternative = "g")
      #p_value_gene[j]=test_output$p.value
      return(test_output$p.value)


    }
    else{
      return(1)
    }

    #return((p_value_gene[j]))
    #return(test_output$p.value)
  }

  #p_value_gene=unlist(lapply(seq(1:length(sub_feature)),perform_sub_entropy))
  p_value_gene=unlist(lapply(seq_len(length(sub_feature)),perform_sub_entropy))

  if(sum(p_value_gene<p_value)>=local_region){
    p_value_final=min(p_value_gene)}
  else{
    p_value_final=1
  }

  return((p_value_final))
}


CIARA=function(norm_matrix,knn_matrix,background,cores_number,p_value,odds_ratio,local_region,approximation){
  run_loop_genes=function(i){
    gene_expression=as.vector(norm_matrix[background[i],])
    print(paste("Gene:",background[i],"number:",i,sep=" "))
    return(CIARA_gene(norm_matrix,knn_matrix,gene_expression,p_value,odds_ratio,local_region,approximation))

  }



 # result_final=do.call(rbind, mclapply(seq(1:length(background)),run_loop_genes , mc.cores = cores_number))

  result_final=do.call(rbind, mclapply(seq_len(length(background)),run_loop_genes , mc.cores = cores_number))


  row.names(result_final)=background

  return(result_final)
}


load(paste0(opt$'norm_matrix',".Rda"))



load(paste0(opt$'knn_matrix',".Rda"))


background=Get_background_full(get(opt$'norm_matrix'),opt$'threshold',opt$'n_cells_low',opt$'n_cells_high')

start_time <- Sys.time()
result=CIARA(get(opt$'norm_matrix'),get(opt$'knn_matrix'),background,opt$'cores_number',opt$'min_p_value',opt$'odds_ratio',opt$'local_region',opt$'approximation')
end_time <- Sys.time()

print(end_time - start_time)





save(result,file=opt$output)



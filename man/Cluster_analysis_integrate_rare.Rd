% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clustering_functions.R
\name{cluster_analysis_integrate_rare}
\alias{cluster_analysis_integrate_rare}
\title{cluster_analysis_integrate_rare}
\usage{
cluster_analysis_integrate_rare(
  raw_counts,
  project_name,
  resolution,
  neighbors,
  max_dimension,
  feature_genes = NULL
)
}
\arguments{
\item{raw_counts}{Raw count matrix (n_genes X n_cells).}

\item{project_name}{Character name of the Seurat project.}

\item{resolution}{Numeric value specifying the parameter \emph{resolution}
used in the Seurat function \emph{FindClusters}.}

\item{neighbors}{Numeric value specifying the parameter \emph{k.param} in
the Seurat function \emph{FindNeighbors}}

\item{max_dimension}{Numeric value specifying the maximum number of the PCA
dimensions used in the parameter \emph{dims} for the Seurat function
\emph{FindNeighbors}}

\item{feature_genes}{vector of features specifying the argument
\emph{features} in the Seurat function \emph{RunPCA}.}
}
\value{
Seurat object including raw and normalized counts matrices, UMAP coordinates and cluster result.
}
\description{
cluster_analysis_integrate_rare
}
\seealso{
\url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindClusters}
\url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindNeighbors}
\url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/RunPCA}
}
\author{
Gabriele Lubatti \email{gabriele.lubatti@helmholtz-muenchen.de}
}

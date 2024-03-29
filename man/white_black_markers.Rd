% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clustering_functions.R
\name{white_black_markers}
\alias{white_black_markers}
\title{white_black_markers}
\usage{
white_black_markers(
  cluster,
  single_cluster,
  norm_matrix,
  marker_list,
  threshold = 0
)
}
\arguments{
\item{cluster}{Vector of length equal to the number of cells, with cluster
assignment.}

\item{single_cluster}{Character. Label of one specify cluster}

\item{norm_matrix}{Norm count matrix (n_genes X n_cells).}

\item{marker_list}{Third element of the output list as returned by the
function \emph{markers_cluster_seurat}}

\item{threshold}{Numeric. The median of the genes across cells belong to
\emph{single_cluster} has to be greater than \emph{threshold} in order to be
consider as a white-black marker for \emph{single_cluster}}
}
\value{
Logical vector of length equal to \emph{marker_list}, with
TRUE/FALSE if the gene is/is not a white-black marker for
\emph{single_cluster}.
}
\description{
A white-marker is a gene whose median expression across cells belong to
\emph{single_cluster} is greater than \emph{threshold} and in all the other
clusters is equal to zero.
}
\author{
Gabriele Lubatti \email{gabriele.lubatti@helmholtz-muenchen.de}
}

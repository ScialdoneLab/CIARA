% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIARA_core_functions.R
\name{get_background_full}
\alias{get_background_full}
\title{get_background_full}
\usage{
get_background_full(
  norm_matrix,
  threshold = 1,
  n_cells_low = 3,
  n_cells_high = 20
)
}
\arguments{
\item{norm_matrix}{Norm count matrix (n_genes X n_cells).}

\item{threshold}{threshold in expression for a given gene}

\item{n_cells_low}{minimum number of cells where a gene is expressed at a
level above threshold}

\item{n_cells_high}{maximum number of cells where a gene is expressed at a
level above threshold}
}
\value{
Character vector with all genes expressed at a level higher than
\emph{threshold} in a number of cells between \emph{n_cells} and
\emph{n_cells_high}.
}
\description{
get_background_full
}
\author{
Gabriele Lubatti \email{gabriele.lubatti@helmholtz-muenchen.de}
}

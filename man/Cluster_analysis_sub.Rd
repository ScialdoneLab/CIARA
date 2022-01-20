\name{Cluster_analysis_sub}
\alias{Cluster_analysis_sub}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cluster_analysis_sub
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Cluster_analysis_sub(raw_counts,resolution,neighbors,max_dimension,name_cluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{raw_counts}{
Raw count matrix (n_genes X n_cells).
}

\item{resolution}{
Numeric value specifying the paramter \emph{resolution} used in the Seurat function \emph{FindClusters}.
}
\item{neighbors}{
Numeric value specifying the parameter \emph{k.param} in the Seurat function \emph{FindNeighbors}
}
\item{max_dimension}{
Numeric value specifying the maximum number of the PCA dimensions used in the  parameter \emph{dims} for the Seurat function \emph{FindNeighbors}
}
\item{name_cluster}{
Character.Name of the original cluster for which the sub clustering is done.

}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Seurat object
}
\references{
\link{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/RunPCA}
\link{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindVariableFeatures}
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

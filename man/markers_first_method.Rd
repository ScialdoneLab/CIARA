\name{markers_first_method}
\alias{markers_first_method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
markers_first_method
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
markers_first_method(seurat_object,clusters,cell_names,number_top)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{seurat_object}{
Seurat object as returned by \emph{Cluster_analysis_sub} or by \emph{Cluster_analysis_integrate_rare}.
}
\item{clusters}{
Vector or length equal to the number of cells, with cluster assignment.
}

\item{cell_names}{
Vector or length equal to the number of cells, with cell names.
}
\item{number_top}{
Integer. Number of top marker genes to keep for each cluster.
}


}
\details{
The Seurat function \emph{FindMarkers} is used to identify general marker for each cluster (specific cluster vs all other cluster). This list of markers is then filtered keeping only the genes that appear as markers in a unique cluster.
}
\value{
List of three elements. The first is a vector with \emph{number_top} marker genes for each cluster. The second is a vector with \emph{number_top} marker genes and corresponding cluster. The third element is a vector with all marker genes for each cluster
}
\references{
%% ~put references to the literature/web site here ~
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
\link{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindMarkers}
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

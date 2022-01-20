\name{Cluster_analysis_integrate_rare}
\alias{Cluster_analysis_integrate_rare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cluster_analysis_integrate_rare
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Cluster_analysis_integrate_rare(raw_counts,batch_name_1,resolution,neighbors,max_dimension,feature_genes=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{raw_counts}{
Raw count matrix (n_genes X n_cells).
}
\item{batch_name_1}{
Character name of the Seurat project.
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
\item{feature_gene}{
vector of features specifying the argument \emph{features} in the Seurat function  \emph{RunPCA}.

}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Seurat object
}
\references{
\link{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindClusters}
\link{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindNeighbors}
\link{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/RunPCA}

}
\author{
Gabriele Lubatti <gabriele.lubatti@helmholtz-muenchen.de>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

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

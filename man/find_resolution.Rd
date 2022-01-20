\name{find_resolution}
\alias{find_resolution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
find_resolution
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find_resolution(seurat_object,resolution_vector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{seurat_object}{
Seurat object as returned by \emph{Cluster_analysis_integrate_rare}
}
\item{resolution_vector}{
vector with all values of resolution for which the Seurat function \emph{FindClusters} is run
}






}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Clustree object
}
\references{
\link{https://cran.r-project.org/web/packages/clustree/clustree.pdf}
}
\author{
Gabriele Lubatti <gabriele.lubatti@helmholtz-muenchen.de>
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

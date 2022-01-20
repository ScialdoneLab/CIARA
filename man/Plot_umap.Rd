\name{Plot_umap}
\alias{Plot_umap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot_umap
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Plot_umap(cordinate_umap,cluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{cordinate_umap}{
Data frame with dimensionality reduction coordinates. Number of rows must be equal to the number of cells
}


\item{cluster}{
Vector with cluster assignment of length equal to the number of cells.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot2 object.
}
\references{
\link{https://cran.r-project.org/web/packages/ggplot2/index.html}
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

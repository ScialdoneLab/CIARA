\name{Plot_genes_sum}
\alias{Plot_genes_sum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot_genes_sum
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Plot_genes_sum(cordinate_umap,norm_counts,genes_relevant,name_title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{cordinate_umap}{
Data frame with dimensionality reduction coordinates. Number of rows must be equal to the number of cells
}
\item{norm_counts}{
Norm count matrix (genes X cells).
}

\item{genes_relevant}{
Vector with gene names for which we want to visualize the sum in each cell.
}
\item{name_title}{
Character value.
}


}
\details{
The sum of each gene in \emph{genes_relevant} across all cells is first normalized to 1.
Then for each cell, the sum from the (normalized) genes expression is computed and shown in the output plot.
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

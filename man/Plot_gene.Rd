\name{Plot_gene}
\alias{Plot_gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot_gene
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Plot_gene(norm_counts,cordinate_umap,gene_id,title_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{norm_counts}{
Norm count matrix (n_genes X n_cells).
}
\item{cordinate_umap}{
Data frame with dimensionality reduction coordinates. Number of rows must be equal to the number of cells
}

\item{gene_id}{
Character name of the gene.
}
\item{title_name}{
Character name.
}


}
\details{
Cells are coloured according to the expression of \emph{gene_id} and plotted according to \emph{cordinate_umap}.
}
\value{
ggplot2 object
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

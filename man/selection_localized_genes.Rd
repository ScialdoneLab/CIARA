\name{selection_localized_genes}
\alias{selection_localized_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
selection_localized_genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
selection_localized_genes(norm_counts,relevant_genes,min_number_cells=4,max_number_genes=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{norm_counts}{
Norm count matrix (n_genes X n_cells).
}
\item{relevant_genes}{
vector of higly localized genes as provided by the last element of the list given as output from \emph{Entropy_mixing_final}.
}
\item{min_number_cells}{
Minimum number of cells where a genes must be expressed (> 0).
}
\item{max_number_genes}{
Maximum number of genes to show for each cell in the interactive plot from \emph{interactive_plot}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Character vector where each entry contains the name of the top \emph{max_number_genes} for the corresponding cell.
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

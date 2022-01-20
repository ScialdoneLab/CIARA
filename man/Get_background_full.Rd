\name{Get_background_full}
\alias{Get_background_full}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get_background_full
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Get_background_full(norm_matrix,treshold,n_cells_low,n_cells_high)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{norm_matrix}{
Norm count matrix (n_genes X n_cells).
}
\item{treshold}{
threshold in expression for a given gene
}

\item{n_cells_low}{
minimum number of cells where a gene is expressed at a level above threshold


}
\item{n_cells_high}{
maximum number of cells where a gene is expressed at a level above threshold
}



}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
vector with all genes expressed at a level higher than \emph{treshold} in a number of cells between \emph{n_cells} and \emph{n_cells_high}
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

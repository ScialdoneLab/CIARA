\name{interactive_plot}
\alias{interactive_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
interactive_plot
}
\description{
It shows in an interactive plot which are the highly localized genes in each cell. It is based on plotly library
}
\usage{
interactive_plot(coordinate_umap,Color,text,min_x=NULL,max_x=NULL,min_y=NULL,max_y=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinate_umap}{
Matrix with dimensionality reduction coordinates, with cells as rows.
}
\item{Color}{
vector of length equal to n_rows in coordinate_umap.Each cell will be coloured following a gradient according to the corresponding value of this vector.
}
\item{text}{
Charachter vector specifying the highly localized genes in each cell. It is the output from \emph{selection_localized_genes}.
}
\item{min_x}{
Set the min limit on the x axis.
}
\item{max_x}{
Set the max limit on the x axis.

}
\item{min_y}{
Set the min limit on the y axis.
}

\item{max_y}{
Set the min limit on the y axis.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
\link{https://plotly.com/r/}

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

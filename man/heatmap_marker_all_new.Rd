\name{heatmap_marker_all_new}
\alias{heatmap_marker_all_new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
heatmap_marker_all_new
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
heatmap_marker_all_new(marker_plot,marker_plot_plot,my.clusters,Condition,norm_counts,number)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{marker_plot}{
First element returned by \emph{markers_first_method}
}
\item{marker_plot_plot}{
Second element returned by \emph{markers_first_method}
}

\item{my.clusters}{
Vector or length equal to the number of cells, with cluster assignment.
}
\item{Condition}{
Vector or length equal to the number of cells, specifying the condition of the cells (i.e. batch, dataset of origin..)
}

\item{norm_counts}{
Norm count matrix (genes X cells).
}

\item{number}{
Size of the text in the heatmap plot.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A Heatmap-class object.
}
\references{
\link{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
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

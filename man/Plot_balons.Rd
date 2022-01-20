\name{Plot_balons}
\alias{Plot_balons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot_balons
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Plot_balons(norm_counts_ane_24,final_cluster_version_sub,marker_plot_day_24_first_plot,max_number,total_cells,max_size=5,text_size=7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{norm_counts_ane_24}{
Normalized count matrix with genes on the rows and cells on the columns.
}
\item{final_cluster_version_sub}{
Vector with cluster information
}
\item{marker_plot_day_24_first_plot}{
Third element of the output list as returned by the function \emph{markers_first_method}
}

\item{max_number}{
Integer. Maximum number of markers for each cluster for which we want to plot the expression.

}

\item{max_size}{
Integer. Size of the dots to be plotted.
}

\item{text_size}{
Integer. Size of the text (legend, x axis and y axis).
}



}
\details{

}
\value{
ggplot object.
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

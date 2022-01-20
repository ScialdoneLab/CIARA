\name{white_black_marker}
\alias{white_black_marker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
white_black_marker
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
white_black_marker(cluster_race_entropy,single_cluster,norm_race,marker_list,threshold=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster_race_entropy}{
Vector with cluster information
}
\item{single_cluster}{
Character. Label of one specifi cluser
}
\item{norm_race}{
Normalized count matrix with genes on the rows and cells on the columns.
}

\item{marker_list}{
Third element of the output list as returned by the function \emph{markers_first_method}

}

\item{threshold}{
Numeric. The median of the genes across cells belong to \emph{single_cluster}  has to be greater than \emph{threshold} in order to be consider a white-black marker for \emph{single_cluster}
}




}
\details{
A white-marker is a gene whose median expression across cells belong to \emph{single_cluster} is greater than \emph{threshold} and in all the other clusters is equal to zero.

}
\value{
Logical vector of length equal to \emph{marker_list}, with TRUE/FALSE if the gene is/is not a white-black marker for \emph{single_cluster}.
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

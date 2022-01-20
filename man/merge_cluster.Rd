\name{merge_cluster}
\alias{merge_cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
merge_cluster
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
merge_cluster(old_cluster,new_cluster,max_number=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{old_cluster}{
original cluster asignment that need to be updated
}
\item{new_cluster}{
new cluster asignment that need to be integrated with \emph{old_cluster}.
}
\item{max_number}{
Threshold in size for clusters in \emph{new_cluster}. Only cluster with number of cells smaller than \emph{max_number} will be integrated in \emph{old cluster}.
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

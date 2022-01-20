\name{test_hvg}
\alias{test_hvg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
test_hvg
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
test_hvg(raw_counts,final_cluster_version,relevant_genes,background,number_hvg,min_p_value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{raw_counts}{
Raw count matrix (n_genes X n_cells).
}
\item{final_cluster_version}{
Vector with cluster assignment for each cell.
}
\item{relevant_genes}{
Character vector with gene names.
}
\item{background}{
Character vector with all the genes names. It is used as backgroud for the Fisher test.
}
\item{number_hvg}{
Integer value. Number of top HVGs provided by the Seurat function \emph{FindVariableFeatures}.
}
\item{min_p_value}{
Threshold on p values provided by Fisher test.
}
}
\details{
For each cluster in \emph{final_cluster_version}, HVGs are defined with Seurat function \emph{FindVariableFeatures}.
A Fisher test is performed to see if there is a statistically significant enrichment between the top \emph{number_hvg} and the \emph{relevant_genes}
}
\value{
A list with two elements.
\item{first element}{The first one is a list with length equal to the number of clusters. Each entry is list of three elements.
The first  two elements contain the p value and the odds ration given by the Fisher test
The third is a vector with genes names that are present both in \emph{relevant_genes} and in top \emph{number_hvg} HVGs
}
\item{second element}{a character vector with the name of the cluster that have a p value smaller than \emph{min_p_value}}
}
\references{
\link{https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/fisher.test}
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

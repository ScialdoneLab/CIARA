\name{Entropy_mixing}
\alias{Entropy_mixing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Entropy_mixing
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Entropy_mixing(norm_matrix,knn_matrix,background,cores_number=1,p_value=0.001,odds_ratio=2,local_region=1,approximation=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{norm_matrix}{
Norm count matrix (n_genes X n_cells).
}
\item{knn_matrix}{
K-nearest neighbors matrix (n_cells X n_cells).
}

\item{background}{
Vector of genes form which the entropy of mixing is computed
}

\item{cores_number}{
Integer.Number of cores to use.
}
\item{p_value}{
p value returned by the function \emph{fisher.test} with parameter alternative = "g"

}

\item{odds_ratio}{
odds_ratio returned by the function \emph{fisher.test} with parameter alternative = "g"

}

\item{local_region}{
Integer. Minimum numer of local regions (cell with its knn neighbours) where the binarized gene expression is enriched in 1.

}
\item{approximation}{
Logical.For a given gene, compute entropy in the local regions of only the cells where the binarized gene expression is 1.
}



}
\details{
For each gene of \emph{background}, it compute entropy of mixing as specified in \emph{Entropy_mixing_gene}
}
\value{
Return a dataframe with n_rows equal to the length of \emph{background} . Each row is the output from \emph{Entropy_mixing_gene}
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

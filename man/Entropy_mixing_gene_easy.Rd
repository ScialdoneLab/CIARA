\name{Entropy_mixing_gene}
\alias{Entropy_mixing_gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Entropy_mixing_gene
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Entropy_mixing_gene(norm_matrix,knn_matrix,gene_expression,p_value=0.001,odds_ratio=2,local_region=1,approximation=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{norm_matrix}{
Norm count matrix (n_genes X n_cells).
}
\item{knn_matrix}{
K-nearest neighbors matrix (n_cells X n_cells).
}
\item{gene_expression}{
numeric vector with the gene expression (length equal to n_cells).
For computing entropy, the gene expression is binarized (equal to 0/1 in the cells where the value is below/above the median)
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
The gene expression is binarized (1/0) if the value in a given cell is above/below the median.
Each of cell with its first K nearest neighbors defined a local region.
If there are at least \emph{local_region} enriched in 1 according to \emph{fisher.test}, then the entropy for the gene is computed starting from the probabilty of having 1/0.
The minimum of the entropy across all the enriched local regions is the entropy of mixing.
If there are no enriched local regions, then the entropy of mixing  and the p value by default are set to 1
}
\value{
List with two elements: the first is the entropy of mixing and the second is the p value
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
 See Also {\link{https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/fisher.test}}

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

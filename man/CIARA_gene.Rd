\name{CIARA_gene}
\alias{CIARA_gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CIARA_gene
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CIARA_gene(norm_matrix,knn_matrix,gene_expression,p_value=0.001,odds_ratio=2,local_region=1,approximation=FALSE)
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
The gene expression is binarized (equal to 0/1 in the cells where the value is below/above the median)
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
Logical.For a given gene, the fisher test is run in the local regions of only the cells where the binarized gene expression is 1.
}




}
\details{
The gene expression is binarized (1/0) if the value in a given cell is above/below the median.
Each of cell with its first K nearest neighbors defined a local region.
If there are at least \emph{local_region} enriched in 1 according to \emph{fisher.test}, then the gene is defined as highly localized and a final p value is assigned to it.
The final p value is the minimum of the p values from all the enriched local regions.
If there are no enriched local regions, then the p value by default is set to 1
}
\value{
List with one element corresponding to the p value of the gene.
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

#' QUBIC: A Qualitative Biclustering Algorithm for Analyses of Gene Expression Data
#'
#' \code{BCQU} performs a QUalitative BIClustering.
#' 
#' For a given representing matrix of a microarray data set, 
#' we construct a weighted graph G with genes represented as vertices, edges connecting every pair of genes, 
#' and the weight of each edge being the similarity level between the two corresponding (entire) rows. 
#' Clearly, the higher a weight, the more similar two corresponding rows are. 
#' Intuitively, genes in a bicluster should induce a heavier subgraph of G because under a subset of the conditions, 
#' these genes have highly similar expression patterns that should make the weight of each involved edge heavier, 
#' comparing to the edges in the background. 
#' But it should be noted that some heavy subgraph may not necessarily correspond to a bicluster, 
#' i.e. genes from a heavy subgraph may not necessarily have similar expression patterns 
#' because different edges in a subgraph may have heavier weights under completely different subsets of conditions. 
#' It should also be noted that recognizing all heavy subgraphs in a weighted graph itself is 
#' computationally intractable because identification of maximum cliques in a graph is a special case of this, 
#' and the maximum clique problem is a well known intractable problem (NP-hard). 
#' So in our solution, we do not directly solve the problem of finding heavy subgraphs in a graph. 
#' Instead, we built our biclustering algorithm based on this graph representation of a microarray gene expression data, 
#' and tackle the biclustering problem as follows. 
#' We find all feasible biclusters (I,J) in the given data set such that min{|I|, |J|} is as large as possible, 
#' where I and J are subsets of genes and conditions, respectively.
#' @name BCQU
#' 
#' @aliases QUBIC QU BCQU-class biclust,matrix,BCQU-method
#'
#' @usage ## S4 method for class 'matrix,BCQU':
#' biclust(x, method = BCQU(), r = 1, q = 0.06, c = 0.95, o = 100, f = 1, k = 2)
#' 
#' @param x The data matrix where biclusters have to be found. 
#' (for example: a qualitative representation of gene expression data)
#' @param r The range of possible ranks.
#' @param q The percentage of the regulating conditions for each gene.
#' @param c The required consistency level of a bicluster.
#' @param o The number of output biclusters.
#' @param f The filter cut-off for data post-processing.
#' @param k The minimum column width of the block, 
#' default: 5\% of columns, minimum 2 columns.
#' 
#' @return Returns an Biclust object, which contains bicluster candidates
#' 
#' @keywords qubic biclust bicluster bi-cluster biclustering bi-clustering
#'
#' @examples
#' #Random matrix with embedded bicluster
#' test <- matrix(rnorm(5000),100,50)
#' test[11:20,11:20] <- rnorm(100,3,0.3)
#' res<-biclust(test, method = BCQU(), r = 1, q = 0.06, c = 0.95, o = 100, f = 1, k = 2)
#' res
#' 
#' \dontrun{
#' #microarray matrix
#' data(BicatYeast)
#' res<-biclust(BicatYeast, method=BCQU())
#' res} 
#' 
#' @references Li G, Ma Q, Tang H, Paterson AH, Xu Y. 
#' QUBIC: a qualitative biclustering algorithm for analyses of gene expression data. 
#' \emph{Nucleic Acids Research}. 2009;\bold{37(15)}:e101. doi:10.1093/nar/gkp491.
NULL

qubiclust <- function(x, r = 1, q = 0.06, c = 0.95, o = 100, f = 1) {
  MYCALL <- match.call()
  res <- qubic(x, r, q, c, o, f)
  return(BiclustResult(as.list(MYCALL),
                       matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
                       matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE),
                       as.numeric(res["Number"]),
                       res["info"]))
}
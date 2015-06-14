#' QUBIC: A Qualitative Biclustering Algorithm for Analyses of Gene Expression Data
#'
#' QUBIC is a biclusting package.
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
#' @name QUBIC
#' 
#' @aliases QUBIC qubic BCQU bcqu BCQU.d bcqu.d biclust method
#'
#' @param x the data matrix where biclusters have to be found. 
#' (for example: a qualitative representation of gene expression data) \cr
#' For \code{BCQU()}, the data matrix should be real \cr
#' For \code{BCQU.d()}, the data matrix should be discretized as integer. 
#' Zeros in the matrix will be treated as non-relevant value.
#' @param r \code{BCQU()} only. The range of possible ranks.
#' @param q \code{BCQU()} only. The percentage of the regulating conditions for each gene.
#' @param c The required consistency level of a bicluster.
#' @param o The number of output biclusters.
#' @param f The filter cut-off for data post-processing.
#' @param k The minimum column width of the block, minimum \code{2} columns.
#' @param P The flag to enlarge current biclsuter by the \emph{p} value constrain.
#' @param S The flag using area as the value of bicluster to determine when stop.
#' @param C The flag using the lower bound of condition number (5 persents of the gene number).
#' 
#' @return Returns an Biclust object, which contains bicluster candidates
#' 
#' @seealso \code{\link{biclust}}
#' 
#' @references Li G, Ma Q, Tang H, Paterson AH, Xu Y. 
#' QUBIC: a qualitative biclustering algorithm for analyses of gene expression data. 
#' \emph{Nucleic Acids Research}. 2009;\bold{37(15)}:e101. doi:10.1093/nar/gkp491.
#' 
#' @keywords qubic biclust bicluster bi-cluster biclustering bi-clustering
NULL

#' \code{BCQU} performs a QUalitative BIClustering.
#' 
#' @name BCQU
#'  
#' @rdname QUBIC
#'
#' @examples
#' #Random matrix with embedded bicluster
#' test <- matrix(rnorm(5000),100,50)
#' test[11:20,11:20] <- rnorm(100,3,0.3)
#' res<-biclust(test, method = BCQU(), r = 1, q = 0.06, c = 0.95, o = 100, f = 1, k = 2, 
#'              P = FALSE, S = FALSE, C = FALSE)
#' res
#'  
#' \dontrun{
#' #Bicluster on microarray matrix
#' data(BicatYeast)
#' res<-biclust(BicatYeast, method=BCQU())
#' res}
#' 
#' \dontrun{
#' #Bicluster on selected of genes
#' data(EisenYeast)
#' res<-biclust(EisenYeast[c("YHR051W","YGL117W","YDR495C"),], method=BCQU())
#' res}
#' 
setClass('BCQU',
         contains = 'BiclustMethod',
         prototype = prototype(
           biclustFunction = function(x,...){.qubiclust(x,...)}))

#' @describeIn QUBIC Performs a QUalitative BIClustering.
#' @usage ## S4 method for class 'matrix,BCQU':
#' biclust(x, method = BCQU(), r = 1, q = 0.06, c = 0.95, o = 100, f = 1, k = 2, 
#'         P = FALSE, S = FALSE, C = FALSE)
BCQU <- function() {
  return(new('BCQU'))
}

#' QUBICD
#' 
#' \code{BCQU.d} performs a QUalitative BIClustering for discretized matrix.
#' 
#' @name BCQU.d-class
#'  
#' @aliases QUBICD QUD BCQU.d-class biclust,matrix,BCQU.d-method
#' 
#' @rdname QUBIC
#'
#' @examples
#' 
#' #
#' #Discretize yeast microarray data
#' data(BicatYeast)
#' res <- biclust(discretize(BicatYeast[1:10,1:10]), method=BCQU.d())
#' res
#' 
setClass('BCQU.d',
         contains = 'BiclustMethod',
         prototype = prototype(
           biclustFunction = function(x,...){.qubiclust_d(x,...)}))

#' @describeIn QUBIC Performs a QUalitative BIClustering for discretized matrix.
#' 
#' @usage ## S4 method for class 'matrix,BCQU.d':
#' biclust(x, method = BCQU.d(), c = 0.95, o = 100, f = 1, k = 2, 
#'         P = FALSE, S = FALSE, C = FALSE)
BCQU.d <- function() {
  return(new('BCQU.d'))
}

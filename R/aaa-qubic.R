######################################################################
#' QUBIC: A Qualitative Biclustering Algorithm for Analyses of Gene Expression Data
#'
#' QUBIC is a biclusting package.
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
#'
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
#' @param type The constrain type.
#' @param P The flag enlarge current biclsuter by the \emph{p} value constrain.
#' @param C The flag using the lower bound of condition number (5 persents of the gene number).
#' @param verbose If '\code{TRUE}', prints extra information on progress.
#'
#' @return Returns an Biclust object, which contains bicluster candidates
#'
#' @details If \code{type} is omitted or \code{type="default"}, the default method of QUBIC will be used.
#'
#' If \code{type="area"}, using area as the value of bicluster to determine when stop.
#'
#' Other types are reserved for future use.
#'
#' @seealso \code{\link{qudiscretize}} \code{\link{qnetwork}} \code{\link{qnet2xml}} \code{\link{biclust}}
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
#' res<-biclust(test, method = BCQU())
#' summary(res)
#' show(res)
#' names(attributes(res))
#'
#' \dontrun{
#' #Load microarray matrix
#' data(BicatYeast)
#' 
#' #Display number of column and row of BicatYeast
#' ncol(BicatYeast)
#' nrow(BicatYeast) 
#' #Bicluster on microarray matrix
#' system.time(res<-biclust(BicatYeast, method=BCQU()))
#' 
#' #Show bicluster info
#' res
#' #Show the first bicluster
#' bicluster(BicatYeast, res, 1)
#' #Get the 4th bicluster
#' bic4 <- bicluster(BicatYeast, res, 4)[[1]]
#' 
#' #or
#' bic4 <- bicluster(BicatYeast, res)[[4]]
#' #Show rownames of the 4th bicluster
#' rownames(bic4)
#' #Show colnames of the 4th bicluster
#' colnames(bic4)
#'
#' }
#' \dontrun{
#' #Bicluster on selected of genes
#' data(EisenYeast)
#' genes <- c("YHR051W", "YKL181W", "YHR124W", "YHL020C", "YGR072W",
#'            "YGR145W", "YGR218W", "YGL041C", "YOR202W", "YCR005C")
#' #same result as res<-biclust(EisenYeast[1:10,], method=BCQU())
#' res<-biclust(EisenYeast[genes,], method=BCQU())
#' res
#'
#' }
#' \dontrun{
#' # Get bicluster by row name = 249364_at
#' bicluster(BicatYeast, res, which(res@@RowxNumber[which(rownames(BicatYeast)=="249364_at"),]))
#'
#' }
#' \dontrun{
#' # Get bicluster by col name = cold_roots_6h
#' bicluster(BicatYeast, res, which(res@@NumberxCol[,which(colnames(BicatYeast)=="cold_roots_6h")]))
#'
#' }
#' \dontrun{
#' #
#' bicluster(BicatYeast, res, which(res@@NumberxCol[,which(colnames(BicatYeast)=="cold_roots_6h")]))
#'
#' }
#' \dontrun{
#' # draw a single bicluster using drawHeatmap {bicust}
#' data(BicatYeast)
#' res <- biclust(BicatYeast, BCQU(), verbose = FALSE)
#' # Draw heatmap of the first cluster
#' drawHeatmap(BicatYeast, res, 1)
#'
#' }
#' \dontrun{
#' # draw a single bicluster using heatmap {stats}
#' data(BicatYeast)
#' res <- biclust(BicatYeast, BCQU(), verbose = FALSE)
#' bic10 <- bicluster(BicatYeast, res, 10)[[1]]
#' 
#' # Draw heatmap of the 10th cluster using heatmap {stats}
#' heatmap(as.matrix(t(bic10)), Rowv = NA, Colv = NA, scale = "none")
#' 
#' # Draw heatmap of the 10th cluster using plot_heatmap {phyloseq}
#' stopifnot(require("phyloseq"))
#' plot_heatmap(otu_table(bic10, taxa_are_rows = TRUE))
#'
#' }
#' \dontrun{
#' # draw a single bicluster with original data background and color options
#' data(BicatYeast)
#' res <- biclust(BicatYeast, BCQU(), verbose = FALSE)
#' palette <- colorRampPalette(c("red", "yellow", "green"))(n = 100)
#' # Draw heatmap of the first cluster with color
#' drawHeatmap(BicatYeast, res, 1, FALSE, beamercolor = TRUE, paleta = palette)
#'
#' }
#' \dontrun{
#' # draw some overlapped biclusters
#' data(BicatYeast)
#' res <- biclust(BicatYeast, BCQU(), verbose = FALSE)
#' biclusternumber(res, 1)
#' biclusternumber(res, 3)
#' # Draw overlapping heatmap
#' heatmapBC(x = BicatYeast, bicResult = res, number = c(1, 3), local = TRUE)
#'
#' }
#' \dontrun{
#' # draw all the biclusters
#' data(BicatYeast)
#' res <- biclust(BicatYeast, BCQU(), verbose = FALSE)
#' # Draw the first bicluster on heatmap
#' heatmapBC(x = BicatYeast, bicResult = res, number = 1)
#' # Draw all the biclusters, not working well.
#' # Overlap plotting only works for two neighbor bicluster defined by the order in the number slot.
#' heatmapBC(x = BicatYeast, bicResult = res, number = 0)
#'
#' }
setClass('BCQU',
         contains = 'BiclustMethod',
         prototype = prototype(
           biclustFunction = function(x,...) {
             .qubiclust(x,...)
           }
         ))

#' @describeIn QUBIC Performs a QUalitative BIClustering.
#' @usage ## S4 method for class 'matrix,BCQU':
#' biclust(x, method = BCQU(), r = 1, q = 0.06, c = 0.95, o = 100, f = 1, k = 2,
#'         type = "default", P = FALSE, C = FALSE, verbose = TRUE)
BCQU <- function() {
  return(new('BCQU'))
}

#' QUBICD
#'
#' \code{BCQU.d} performs a QUalitative BIClustering for a discret matrix.
#'
#' @name BCQU.d-class
#'
#' @aliases qubic_d QUBICD QUD BCQU.d-class biclust,matrix,BCQU.d-method
#'
#' @rdname QUBIC
#'
#' @examples
#' #Biclustering of discretized yeast microarray data
#' data(BicatYeast)
#' disc<-qudiscretize(BicatYeast[1:10,1:10])
#' biclust(disc, method=BCQU.d())
setClass('BCQU.d',
         contains = 'BiclustMethod',
         prototype = prototype(
           biclustFunction = function(x,...) {
             .qubiclust_d(x,...)
           }
         ))

#' @describeIn QUBIC Performs a QUalitative BIClustering for a discret matrix.
#'
#' @usage ## S4 method for class 'matrix,BCQU.d':
#' biclust(x, method = BCQU.d(), c = 0.95, o = 100, f = 1, k = 2,
#'         type = "default", P = FALSE, C = FALSE, verbose = TRUE)
BCQU.d <- function() {
  return(new('BCQU.d'))
}

#' Create a qualitative discret matrix
#'
#' \code{qudiscretize} delivers a discret matrix.
#'
#' @usage qudiscretize(matrix, r = 1L, q = 0.06)
#' @param r The range of possible ranks.
#' @param q The percentage of the regulating conditions for each gene.
#'
#' @name qudiscretize
#'
#' @aliases qudiscretize
#'
#' @examples
#' #Qualitative discretize yeast microarray data
#' data(BicatYeast)
#' qudiscretize(BicatYeast[1:7, 1:5])
#'
#' @seealso \code{\link{QUBIC}} \code{\link{discretize}} \code{\link{qugraph}}
NULL
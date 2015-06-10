setClass('BCQU',
         contains = 'BiclustMethod',
         prototype = prototype(
           biclustFunction = function(x,...){qubiclust(x,...)}))


#' QUBIC: A Qualitative Biclustering Algorithm for Analyses of Gene Expression Data
#'
#' \code{BCQU} performs a QUalitative BIClustering.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param x The data matrix where biclusters have to be found. 
#' (for example: a qualitative representation of gene expression data)
#' @return Returns an Biclust object.
#' @references Li G, Ma Q, Tang H, Paterson AH, Xu Y. QUBIC: a qualitative biclustering algorithm for analyses of gene expression data. \emph{Nucleic Acids Research}. 2009;\bold{37(15)}:e101. doi:10.1093/nar/gkp491.
#' @examples
#' #Random matrix with embedded bicluster
#' test <- matrix(rnorm(5000),100,50)
#' test[11:20,11:20] <- rnorm(100,3,0.3)
#' res<-biclust(test, method=BCQU())
#' res
#' 
#' \dontrun{
#' #microarray matrix
#' data(BicatYeast)
#' res<-biclust(BicatYeast, method=BCQU())
#' res
#' }
BCQU <- function() {
  return(new('BCQU'))
}

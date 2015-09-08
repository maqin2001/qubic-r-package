#' Create a qualitative discrete matrix for a given gene expression matrix
#'
#' \code{qudiscretize} delivers a discret matrix.
#'
#' @usage qudiscretize(matrix, r = 1L, q = 0.06)
#' @inheritParams QUBIC
#'
#' @name qudiscretize
#'
#' @aliases qudiscretize
#'
#' @examples
#' # Qualitative discretize yeast microarray data
#' data(BicatYeast)
#' qudiscretize(BicatYeast[1:7, 1:5])
#'
#' @seealso \code{\link{QUBIC}} \code{\link{discretize}} \code{\link{qugraph}}
qudiscretize <- function(matrix, r = 1L, q = 0.06) {
  qubic_discretize(matrix, r, q)
}
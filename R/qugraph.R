#' Draw Bicluster Graph
#'
#' This is the main function of \code{qugraph}
#' which automatically creates an appropriate network for one or two biclusters and
#' sends it to the plotting method.
#' @param x The data matrix
#' @param BicRes BiclustResult object
#' @param number Which bicluster to be plotted
#' @param method a character string indicating
#' which correlation coefficient (or covariance) is to be computed.
#' One of "pearson" (default), "kendall", or "spearman", can be abbreviated.
#' @param ... Any additional arguments described in \code{\link{qgraph}}.
#' @examples
#' \dontrun{
#' #Load microarray matrix
#' data(BicatYeast)
#' res<-biclust(BicatYeast[1:50, ], method=BCQU(), verbose = FALSE)
#' #Draw two biclusters
#' qugraph(BicatYeast[1:50, ], res, number = c(4, 13), group = c(4, 13), method = "spearman",
#'         layout = "spring", minimum = 0.6, edge.labe = FALSE)
#' }
#' \dontrun{
#' #Load microarray matrix
#' data(BicatYeast)
#' res<-biclust(BicatYeast[1:50, ], method=BCQU(), verbose = FALSE)
#' #Draw all biclusters
#' qugraph(BicatYeast[1:50, ], res, group = c(4, 13), method = "spearman",
#'         layout = "spring", minimum = 0.6, edge.labe = FALSE)
#' }
#' @seealso \code{\link{QUBIC}} \code{\link{qgraph}} \code{\link{cor}}
qugraph <- function(x, BicRes, number = 1:BicRes@Number,
                    group = c(number[[1]]),
                    method = c("pearson", "kendall", "spearman"),
                    layout = "spring", minimum = 0.6, ...) {
  if (length(number) < 1)
    stop("at least 1 bicluster needed.")
  if (is.null(rownames(x)))
    stop("can not plot without rownames.")
  if (is.null(colnames(x)))
    stop("can not plot without colnames.")
  bics <- bicluster(x, BicRes, number)
  index <- which(number %in% group)
  
  rownamelist <- list()
  colnamelist <- list()
  for (i in 1:length(bics)) {
    rownamelist[[names(bics)[i]]] <- rownames(bics[[i]])
    colnamelist[[names(bics)[i]]] <- colnames(bics[[i]])
  }
  
  allrownames <- Reduce(union, rownamelist)
  allcolnames <- Reduce(union, colnamelist)
  
  un <- x[allrownames, allcolnames]
  rowidlist <- list()
  
  if (length(group) != 2)
    stop("length(group) != 2")
  rowidlist[[paste(names(bics)[index[[1]]], " & ", names(bics)[index[[2]]], sep = "")]] <-
    match(intersect(rownamelist[[index[[1]]]], rownamelist[[index[[2]]]]), rownames(un))
  rowidlist[[names(bics)[index[[1]]]]] <-
    match(setdiff(rownamelist[[index[[1]]]], rownamelist[[index[[2]]]]), rownames(un))
  rowidlist[[names(bics)[index[[2]]]]] <-
    match(setdiff(rownamelist[[index[[2]]]], rownamelist[[index[[1]]]]), rownames(un))
  rowidlist[["Others"]] <-
    match(setdiff(allrownames, union(rownamelist[[index[[1]]]], rownamelist[[index[[2]]]])), rownames(un))
  
  cort <- cor(t(un), method = method)
  
  stopifnot(require("qgraph"))
  qgraph(
    cort, groups = rowidlist, layout = layout, minimum = minimum, color = cbind(rainbow(length(rowidlist) - 1),"gray"), ... = ...
  )
}
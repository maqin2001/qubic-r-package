#' @describeIn QUBIC Performs a QUalitative BIClustering for a discret matrix.
#'
#' @usage qubiclust_d(x, c = 0.95, o = 100, f = 1,
#'         k = max(ncol(x) \%/\% 20, 2),
#'         type = 'default', P = FALSE, C = FALSE, verbose = TRUE,
#'         weight = NULL, seedbicluster = NULL)
qubiclust_d <- function(x, c = 0.95, o = 100, f = 1,
                         k = max(ncol(x)%/%20, 2), type = "default",
                         P = FALSE, C = FALSE, verbose = TRUE, weight = NULL, seedbicluster = NULL) {
  MYCALL <- match.call()
  S <- (type == "area")
  if(!is.null(weight)) {
    w <- Matrix::Matrix(data = 0, nrow = nrow(x), ncol = nrow(x), dimnames = list(rownames(x), rownames(x)))
    weight[] <- rank(weight, ties.method = "average")
    intersect_rows <- intersect(rownames(x), rownames(weight))
    w[intersect_rows, intersect_rows] <- weight[intersect_rows, intersect_rows]
    res <- .qubic_dw(x, c, o, f, k, P, S, C, verbose, w)
  } else if(!is.null(seedbicluster)) {
    if (seedbicluster@Number >= 1)
      res <- .qubic_de(x, c, verbose, seedbicluster@RowxNumber, seedbicluster@NumberxCol)
    else
      return(seedbicluster)
  } else res <- .qubic_d(x, c, o, f, k, P, S, C, verbose)
  return(biclust::BiclustResult(as.list(MYCALL), matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
                                matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
                                res["info"]))
}

#' @describeIn Input scRNA expression data, do biclustering + dense graph, output dense graph.
scqubiclust_d <- function(x, c = 0.95, o = 100, f = 1,
                        k = max(ncol(x)%/%20, 2), type = "default",
                        P = FALSE, C = FALSE, verbose = TRUE, weight = NULL, seedbicluster = NULL) {
  MYCALL <- match.call()
  S <- (type == "area")
  if(!is.null(weight)) {
    w <- Matrix::Matrix(data = 0, nrow = nrow(x), ncol = nrow(x), dimnames = list(rownames(x), rownames(x)))
    weight[] <- rank(weight, ties.method = "average")
    intersect_rows <- intersect(rownames(x), rownames(weight))
    w[intersect_rows, intersect_rows] <- weight[intersect_rows, intersect_rows]
    res <- .scqubic_dw(x, c, o, f, k, P, S, C, verbose, w)
  } else if(!is.null(seedbicluster)) {
    if (seedbicluster@Number >= 1)
      res <- .scqubic_de(x, c, verbose, seedbicluster@RowxNumber, seedbicluster@NumberxCol)
    else
      return(seedbicluster)
  } else res <- .scqubic_d(x, c, o, f, k, P, S, C, verbose)
  rownames(res) <- colnames(x)
  colnames(res) <- colnames(x)
  return(res)
}

#' @describeIn QUBIC Performs a QUalitative BIClustering.
#'
#' @usage qubiclust(x, r = 1L, q = 0.06, c = 0.95, o = 100, f = 1,
#'         k = max(ncol(x) \%/\% 20, 2),
#'         type = 'default', P = FALSE, C = FALSE, verbose = TRUE,
#'         weight = NULL, seedbicluster = NULL)
qubiclust <- function(x, r = 1L, q = 0.06, c = 0.95, o = 100, f = 1,
                       k = max(ncol(x)%/%20, 2), type = "default",
                       P = FALSE, C = FALSE, verbose = TRUE, weight = NULL, seedbicluster = NULL) {
  x_d <- qudiscretize(x, r, q)
  return(qubiclust_d(x_d, c, o, f, k, type, P, C, verbose, weight, seedbicluster))
}

#' @describeIn Input scRNA expression data, do discritization + biclustering + dense graph, output dense graph.
scqubiclust <- function(x, r = 1L, q = 0.06, c = 0.95, o = 100,
  f = 1, k = max(ncol(x)%/%20, 2), type = "default", P = FALSE,
  C = FALSE, verbose = TRUE, weight = NULL, seedbicluster = NULL) {
  x_d <- qudiscretize(x, r, q)
  return(scqubiclust_d(x_d, c, o, f, k, type, P, C, verbose,
    weight, seedbicluster))
}

#' @describeIn Input bicluster, do dense graph, output graph matrix (symmetric)
scgraph <- function(res, cleardiagonal = FALSE, dimnames = NULL) {
  nc <- ncol(res@NumberxCol)
  matrix <- rep(FALSE, nc) %o% rep(FALSE, nc)
  number <- nrow(res@NumberxCol)
  for (index in 1:number) {
    col <- res@NumberxCol[index, ]
    matrix <- matrix + col %o% col
  }
  if (!is.null(dimnames))
    dimnames(matrix) <- dimnames
  if (cleardiagonal)
    diag(matrix) <- 0
  return(matrix)
}

#' @describeIn Input multiple biclusters, do dense graph for each and add them up, output single graph matrix (symmetric)
scmgraph <- function(reslist, cleardiagonal = FALSE, dimnames = NULL) {
  nc <- ncol(reslist[[1]]@NumberxCol)
  matrix <- rep(FALSE, nc) %o% rep(FALSE, nc)
  for (res in reslist) {
    if (nc != ncol(res@NumberxCol))
      stop("Different NumberxCol")
    number <- nrow(res@NumberxCol)
    for (index in 1:number) {
      col <- res@NumberxCol[index, ]
      matrix <- matrix + col %o% col
    }
  }
  if (!is.null(dimnames))
    dimnames(matrix) <- dimnames
  if (cleardiagonal)
    diag(matrix) <- 0
  return(matrix)
}

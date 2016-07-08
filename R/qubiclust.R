.qubiclust_d <- function(x, c = 0.95, o = 100, f = 1,
                         k = max(ncol(x)%/%20, 2), type = "default",
                         P = FALSE, C = FALSE, verbose = TRUE, weight = NULL, ...) {
  MYCALL <- match.call()
  S <- (type == "area")
  if(is.null(weight)) res <- .qubic_d(x, c, o, f, k, P, S, C, verbose, ...)
  else {
    w <- matrix(nrow = ncol(x), ncol = ncol(x), dimnames = list(rownames(x), rownames(x)))
    weight[] <- rank(weight, ties.method = "average")
    for (i in rownames(x))
      for (j in rownames(x))
        if((i %in% rownames(weight)) && (j %in% rownames(weight)))
          w[i,j] <- weight[i,j]
        else
          w[i,j] <- 0
    res <- .qubic_dw(x, c, o, f, k, P, S, C, verbose, w, ...)
  }
  return(biclust::BiclustResult(as.list(MYCALL), matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
                                matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
                                res["info"]))
}

.qubiclust <- function(x, r = 1L, q = 0.06, c = 0.95, o = 100, f = 1,
                       k = max(ncol(x)%/%20, 2), type = "default",
                       P = FALSE, C = FALSE, verbose = TRUE, weight = NULL, ...) {
  x_d <- qudiscretize(x, r, q)
  return(.qubiclust_d(x_d, c, o, f, k, type, P, C, verbose, weight, ...))
}

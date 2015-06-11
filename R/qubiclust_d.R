qubiclust_d <- function(x, c = 0.95, o = 100, f = 1, k = 2) {
  MYCALL <- match.call()
  res <- qubic_d(x, c, o, f, k)
  return(BiclustResult(as.list(MYCALL),
                       matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
                       matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE),
                       as.numeric(res["Number"]),
                       res["info"]))
}
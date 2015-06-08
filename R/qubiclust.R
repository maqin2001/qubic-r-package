qubiclust <- function(x) {
  MYCALL <- match.call()
  res <- qubic(x)
  return(BiclustResult(as.list(MYCALL),
                       matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
                       matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE),
                       as.numeric(res["Number"]),
                       res["info"]))
}
library(QUBIC)

runme <- function(nrow, ncol = nrow) {
  test <- matrix(rnorm(5000), nrow, ncol)
  test[11:20, 11:20] <- rnorm(100, 3, 0.3)
  system.time(res <- biclust::biclust(test, method = BCQU(), verbose = FALSE))
}

runme(100)
#user  system elapsed
#0.01    0.00    0.01
#runme(1000)
#user  system elapsed
#0.42    0.00    0.42
#> runme(2000)
#user  system elapsed
#'3.2     0.0     3.2
#'> runme(3000)
#'user  system elapsed
#'10.54    0.02   10.60
#'> runme(4000)
#'user  system elapsed
#'22.08    0.10   22.25
#'> runme(5000)
#'user  system elapsed
#'37.24    0.08   37.47
#'> runme(6000)

#runme(6000)
#runme(7000)
#runme(8000)
#runme(9000)
#runme(10000)
#runme(15000)
#runme(20000)
#runme(25000)

library(igraph)
file = "testdata/example_weight";
graph = read.graph(file, format = "ncol")
get.edgelist(graph, names=TRUE)
E(graph)$weight

test <- matrix(rnorm(5000), 300, 300)
test[11:20, 11:20] <- rnorm(100, 3, 0.3)
system.time(res <- biclust::biclust(test, method = BCQU(), verbose = FALSE))

showinfo <- function(test, res) {
  cat("number of detected biclusters:")
  cat(res@Number)
  cat("\n")
  bic <- biclust::bicluster(test, res, 1)[[1]]
  bicall <- biclust::bicluster(test, res)
  cat("nrow of the first bicluster:")
  cat(nrow(bic))
  cat("\n")
  cat("ncol of the first bicluster:")
  cat(ncol(bic))
  cat("\n")
  cat("area of the first bicluster:")
  cat(nrow(bic) * ncol(bic))
  cat("\n")
  cat("ratio (nrow / ncol) of the first bicluster:")
  cat(nrow(bic) / ncol(bic))
  cat("\n")
  cat("ratio (nrow / ncol) of the matrix:")
  cat(nrow(test) / ncol(test))
  cat("\n")
  maxnrow <- c(-1, -1)
  maxncol <- c(-1, -1)
  maxarea <- c(-1, -1)
  for (i in 1:res@Number) {
    nrow <- nrow(bicall[[i]])
    ncol <- ncol(bicall[[i]])
    area <- nrow * ncol
    if (maxnrow[[1]] < nrow) maxnrow = c(nrow, i)
    if (maxncol[[1]] < ncol) maxncol = c(ncol, i)
    if (maxarea[[1]] < area) maxarea = c(area, i)
  }
  cat("max nrow and corresponding bicluster:")
  cat(maxnrow)
  cat("\n")
  cat("max ncol and corresponding bicluster:")
  cat(maxncol)
  cat("\n")
  cat("max area and corresponding bicluster:")
  cat(maxarea)
  cat("\n")
  cat("union of rows, (# and %):")
  genes_union <- sum(apply(res@RowxNumber, 1, max))
  cat(c(genes_union, genes_union / nrow(res@RowxNumber) * 100))
  cat("\n")
  cat("union of columns, (# and %):")
  conditions_union <- sum(apply(res@NumberxCol, 2, max))
  cat(c(conditions_union, conditions_union / ncol(res@NumberxCol) * 100))
  cat("\n")
  cat("overlap of first two biclusters (row, col, area):")
  rows_overlap <- sum(apply(res@RowxNumber[,1:2], 1, min))
  cols_overlap <- sum(apply(res@NumberxCol[1:2,], 2, min))
  cat(c(rows_overlap, cols_overlap, rows_overlap * cols_overlap))
  cat("\n")
}

showinfo(test, res)


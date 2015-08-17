#' Convert newwork to XGMML
#' @param net Result of \code{\link{qnetwork}}
#' @param minimum cutoff, default: 0.6
#' @param color default: cbind(rainbow(length(net[[2]]) - 1), "gray")
#' @examples
#' \dontrun{
#' # Load microarray matrix
#' data(BicatYeast)
#' res<-biclust(BicatYeast[1:50, ], method=BCQU(), verbose = FALSE)
#' # Get all biclusters
#' net <- qnetwork(BicatYeast[1:50, ], res, group = c(4, 13), method = "spearman")
#' # Save the network to a XGMML file
#' sink("tempnetworkresult.gr")
#' qnet2xml(net, minimum = 0.6, color = cbind(rainbow(length(net[[2]]) - 1), "gray"))
#' sink()
#' # You can use Cytoscape, Biomax or JNets open file named tempnetworkresult.gr
#' 
#' }
#' @seealso \code{\link{qnetwork}} \code{\link{QUBIC}}
qnet2xml <-
  function(net, minimum = 0.6, color = cbind(rainbow(length(net[[2]]) - 1), "gray")) {
    cat('<?xml version="1.0" encoding="utf-8"?>')
    cat("\n")
    cat('<!DOCTYPE graph SYSTEM "http://www.cs.rpi.edu/~puninj/XGMML/xgmml.dtd">')
    cat("\n")
    cat('<graph directed="0" label="QUGRAPH: Hello, I am a graph" >')
    cat("\n")
    for (i in 1:length(net[[2]])){
      for (j in net[[2]][[i]]) {
        cat('  <node label="node ')
        cat(j)
        cat('" id="')
        cat(j)
        cat('"><graphics type="ELLIPSE" fill="')
        cat(substr(color[[i]],1,7))
        cat('"/></node>')
        cat("\n")
      }
    }
    for (i in 1:(nrow(net[[1]]) - 1)) {
      for (j in (i + 1):nrow(net[[1]])) {
        if (net[[1]][i,j] >= minimum){
          cat('  <edge source="')
          cat(i)
          cat('" target="')
          cat(j)
          cat('" ><att name="weight" type="real" value="')
          cat(net[[1]][i,j])
          cat('"/></edge>')
          cat("\n")
        }
      }
    }
    cat('</graph>')
    cat("\n")
  }
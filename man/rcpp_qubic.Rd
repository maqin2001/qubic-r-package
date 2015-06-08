\name{BCQU}
\alias{BCQU}
\docType{package}
\title{The QU Biclustering algorithm}
\description{The function does greedy biclustering of the data matrix.}
\usage{biclust(x, method=BCQU())}
\arguments{\item{x}{Data matrix.}}
\value{Returns an object of class \code{Biclust}.}
\examples{
test <- matrix(rbinom(400, 50, 0.4), 20, 20)
res <- biclust(test, method=BCQU())
res
}
\keyword{cluster}
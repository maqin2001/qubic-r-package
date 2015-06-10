setClass('BCQU',
         contains = 'BiclustMethod',
         prototype = prototype(
           biclustFunction = function(x,...){qubiclust(x,...)}))

BCQU <- function() {
  return(new('BCQU'))
}

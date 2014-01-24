reads.v.fitness <- function(reads, fitnesses, breaks){
  #create vector to add bins
  bins <- list()
  for( i in seq(1,length(reads)) ){
    ptm <- proc.time()
    for( j in seq(2,length(breaks)) ){
      if( reads[i] <= breaks[j] && reads[i] > breaks[j-1] ){
        bins <- c(bins, j)
      }
    }
    if( i %% 1000 == 0){ print(c(i, proc.time()-ptm)) }
  }
  return(cbind(reads, fitnesses, unlist(bins)))
}
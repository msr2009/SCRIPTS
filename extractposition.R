extract.position <- function(x){
  out <- data.frame(wt="", pos=0, mut="")
  out$wt <- substring(x,1,1)
  out$mut <- substring(x, nchar(x))
  out$pos <- as.numeric(substring(x,2,nchar(x)-1))
  return(out)
  }

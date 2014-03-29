compare.frequencies <- function( file1, file2, knownlist, offset=0, outfile="RPlot" ){
  
  offset <- as.numeric(offset)
  d1 <- read.table(file1, sep="\t", header=F)
  d2 <- read.table(file2, sep="\t", header=F)
    
  known <- read.table(knownlist)
  known <- unlist(known)
  k <- data.frame("pos" <- known)
  names(k) <- "pos"
  
  #calculate frequencies
  d1f <- calc_freq(d1, offset)
  d2f <- calc_freq(d2, offset)
  
  #merge with known list, 
  #replace non-polymorphic sites with 1.0
  #low-coverage sites = NA
  d1f.k <- merge(k, d1f, by="pos", all.x=T)
  d1f.k[is.na(d1f.k)] <- 1.0
  d1f.k[d1f.k == 1000] <- NA
  
  d2f.k <- merge(k, d2f, by="pos", all.x=T)
  d2f.k[is.na(d2f.k)] <- 1.0  
  d2f.k[d2f.k == 1000] <- NA
  
  #plot them
  setEPS()
  postscript(outfile)
  plot(d1f.k$freq, d2f.k$freq, xlab=paste(file1, "ALT frequency", sep=" "), ylab=paste(file2, "ALT frequency", sep=" "), col="blue", xlim=c(0,1), ylim=c(0,1) )
  abline(a=1, b=-1, lty=2, col="red")
#  abline(lm(d1f.k$freq ~ d2f.k$freq), col="red")
  dev.off()
}

calc_freq <- function(bed, offset=0){
  depths <- apply(bed[,5:8], 1, sum)
  f <- (bed$V5+bed$V6)/depths
  d <- data.frame("pos" = bed$V2-offset, "freq" = f)
  d$freq[depths < .1*mean(depths)] <- 1000
  return(d)
}

args <- commandArgs(trailingOnly = T)
compare.frequencies(args[1], args[2], args[3], args[4], args[5])

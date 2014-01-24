#function to plot mutation frequencies

#first get commands
args <- commandArgs(TRUE)
filename = args[1]

#load function to extract position information from annotation
extract.position <- function(x){
  out <- data.frame(wt="", pos=0, mut="")
  out$wt <- substring(x,1,1)
  out$mut <- substring(x, nchar(x))
  out$pos <- as.numeric(substring(x,2,nchar(x)-1))
  return(out)
}

#load stringr library to do fun stuff with strings
library(stringr)

#get filename
f <- strsplit(filename, split='/')[[1]]
f <- tail(f, n=1)
#make sure that the genelength is an integer
l <- as.integer(args[2])
#open BED file
dat <- read.delim(filename, sep='\t')

#make a vector containing the annotation data
dat.pos <- str_split_fixed(dat$annotation, ",", 2)[,1]
dat.annotate <- data.frame(wt=substring(dat.pos,1,1), pos=as.numeric(substring(dat.pos,2,nchar(dat.pos)-1)), mut=substring(dat.pos,nchar(dat.pos)))
#get nonsynonymous, nonsense mutations
nonsynonymous <- which(as.character(dat.annotate$mut) != as.character(dat.annotate$wt))
nonsense <- which(as.character(dat.annotate$mut) == "*")

#open PDF output file
pdf(paste(filename, ".plot_mutations.pdf", sep=""), width=12, height=8)

dat.frames <- dat$start[as.numeric(dat$start) < l & dat$type %in% c("D-1","D-2","D-4","D-5","I-1","I-2","I-4","I-5")]
dat.nons <- dat$start[as.numeric(dat$start) < l & dat.annotate$mut == "*"]
dat.exclude <- c(which(as.numeric(dat$start) < l & dat$type %in% c("D-1","D-2","D-4","D-5","I-1","I-2","I-4","I-5")), which(as.numeric(dat$start) < l & as.character(dat.annotate$mut) == "*"))

#remove excluded points from dataset (for plotting purposes)
dat.subs <- dat$start[-dat.exclude]

plot(tabulate(dat.subs[as.numeric(dat.subs) < l])/length(dat$start), main=f, xlab="Position", ylab="% mutations at position")
points(tabulate(dat.frames)/length(dat$start), col="blue")
points(tabulate(dat.nons)/length(dat$start), col="red")

#plots for protein mutations
#1: all mutations
#2: just nonsynonymous (no nonsense)
#3: just nonsense

dev.off()
  

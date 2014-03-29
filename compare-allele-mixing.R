calc_freq <- function(bed, offset=0){
  depths <- apply(bed[,5:8], 1, sum)
  f <- (bed$V5+bed$V6)/depths
  d <- data.frame("pos" = bed$V2-offset, "freq" = f)
  d$freq[depths < .1*mean(depths)] <- 1000
  return(d)
}

setwd("/Volumes/mattrich/DUNHAM/FLO8/20140318_poolsequencing/")

#known mutations
known <- read.table("/Volumes/mattrich/DUNHAM/FLO8/FLO8-mismatches.txt")
known <- unlist(known)
k <- data.frame("pos" <- known)
names(k) <- "pos"

d1 <- read.table("s_1_1_GGCAATTCT.DBVPG1106+1kb.bed", header=F, sep="\t")
d10 <- read.table("s_1_1_GTTACAGGT.DBVPG1106+1kb.bed", header=F, sep="\t")
d25 <- read.table("s_1_1_ACACTGATG.DBVPG1106+1kb.bed", header=F, sep="\t")
d50 <- read.table("s_1_1_TGAGGATTC.DBVPG1106+1kb.bed", header=F, sep="\t")
d75 <- read.table("s_1_1_GTATACAGC.DBVPG1106+1kb.bed", header=F, sep="\t")
d90 <- read.table("s_1_1_ATAGCCTTC.DBVPG1106+1kb.bed", header=F, sep="\t")
d99 <- read.table("s_1_1_ATCGATTCC.DBVPG1106+1kb.bed", header=F, sep="\t")

expect = c(.01, .1, .25, .5, .75, .9, .99)

d1.f <- calc_freq(d1, 1000)
d10.f <- calc_freq(d10, 1000)
d25.f <- calc_freq(d25, 1000)
d50.f <- calc_freq(d50, 1000)
d75.f <- calc_freq(d75, 1000)
d90.f <- calc_freq(d90, 1000)
d99.f <- calc_freq(d99, 1000)

d.all <- Reduce(function(x,y) merge(x, y, all.x=T, by="pos"), list(k, d1.f, d10.f, d25.f, d50.f, d75.f, d90.f, d99.f))
names(d.all) <- c("pos", "d1.f", "d10.f", "d25.f", "d50.f", "d75.f", "d90.f", "d99.f")
d.all[is.na(d.all)] <- 1.0
d.all[,2:8] <- 1-d.all[,2:8]
d.all
#~~~~~~~~~~

k1 <- read.table("s_1_1_GGCAATTCT.K11+1kb.bed", header=F, sep="\t")
k10 <- read.table("s_1_1_GTTACAGGT.K11+1kb.bed", header=F, sep="\t")
k25 <- read.table("s_1_1_ACACTGATG.K11+1kb.bed", header=F, sep="\t")
k50 <- read.table("s_1_1_TGAGGATTC.K11+1kb.bed", header=F, sep="\t")
k75 <- read.table("s_1_1_GTATACAGC.K11+1kb.bed", header=F, sep="\t")
k90 <- read.table("s_1_1_ATAGCCTTC.K11+1kb.bed", header=F, sep="\t")
k99 <- read.table("s_1_1_ATCGATTCC.K11+1kb.bed", header=F, sep="\t")

k1.f <- calc_freq(k1, 1000)
k10.f <- calc_freq(k10, 1000)
k25.f <- calc_freq(k25, 1000)
k50.f <- calc_freq(k50, 1000)
k75.f <- calc_freq(k75, 1000)
k90.f <- calc_freq(k90, 1000)
k99.f <- calc_freq(k99, 1000)

k.all <- Reduce(function(x,y) merge(x, y, all.x=T, by="pos"), list(k, k1.f, k10.f, k25.f, k50.f, k75.f, k90.f, k99.f))
names(k.all) <- c("pos", "k1.f", "k10.f", "k25.f", "k50.f", "k75.f", "k90.f", "k99.f")
k.all[is.na(k.all)] <- 1.0

d.mean <- apply( d.all[,2:8], 2, mean )
d.sd <- apply( d.all[,2:8] , 2, sd)
k.mean <- apply( k.all[,2:8], 2, mean )
k.sd <- apply( k.all[,2:8] , 2, sd)

plot(expect, expect, xlim=c(0,1), ylim=c(0,1.2), type="n", axes=F, xlab="Expected K11 allele frequecy", ylab="Observed K11 allele frequency")
axis(1, at=c(.01, .1, .25, .5, .75, .9, .99))
axis(2)
arrows(x0=expect, y0=k.mean-k.sd, x1=expect, y1=k.mean+k.sd, angle=90, code=3)
points(expect, k.mean, pch=16, cex=1.5, col="yellow")

arrows(x0=expect, y0=d.mean-d.sd, x1=expect, y1=d.mean+d.sd, angle=90, code=3, col="blue", lty=2)
points(expect, d.mean, pch=16, cex=1.5, col="blue")
abline(h=1.0, lty=3)
legend(x=.85, y=.2, col=c("blue", "yellow"), pch=16, legend=c("DBVGP1106", "K11"))
text(x=.93, y=.22, labels="Alignment to:")

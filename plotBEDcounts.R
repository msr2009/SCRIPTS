args <- commandArgs(TRUE)
infile <- args[1]

d <- read.delim(infile, sep='\t', header=T)

#filter out low read-counts
d <- d[-which(d[,6] <= args[2]),]

pdf(paste(strsplit(infile, split=".txt")[[1]][1], ".plots.pdf", sep=""), width=12, height=8)

#~~Position plot~~~~~~~~~~~~~~~~~~~~#
plot.new()
plot.window(xlim=c(0,max(d$aa_position)), ylim=c(min(d$count/sum(d$count)),max(d$count/sum(d$count))))
axis(1, at=(c(seq(0,max(d$aa_position),by=50), max(d$aa_position))))
axis(2)

#plot synonymous changes
points(d$aa_position[which(d$aa_type=="S")], d$count[which(d$aa_type=="S")]/sum(d$count), col="red", pch=16)
#plot nonsynonymous changes
points(d$aa_position[which(d$aa_type=="NS")], d$count[which(d$aa_type=="NS")]/sum(d$count), col="blue", pch=16)
#plot nonsense changes
points(d$aa_position[which(d$aa_type=="N")], d$count[which(d$aa_type=="N")]/sum(d$count), col="green", pch=16)
#plot frameshifts
points(d$aa_position[which(d$aa_type=="F")], d$count[which(d$aa_type=="F")]/sum(d$count), col="orange", pch=16)

title(main=infile, xlab="AA position", ylab="Percent of total mutations")

dev.off()

#~~Boxplots~~~~~~~~~~~~~~~~~~~~~~#
#boxplot(d$count ~ d$aa_type, col=c("orange", "green", "blue", "red"), main=infile, ylab="log2 Enrichment", xlab="Mutation Type")

#~Scatters~~~~~~~~~~~~~~~~~~~~~~~#
#plot.new()
#plot.window(xlim=c(0,max(d[,7])), ylim=c(0,max(d[,9])))
#axis(1); axis(2)
#points(d[which(d$aa_type=="S"),7], d[which(d$aa_type=="S"),9], col="red", pch=16)
#points(d[which(d$aa_type=="NS"),7], d[which(d$aa_type=="NS"),9], col="blue", pch=16)
#points(d[which(d$aa_type=="N"),7], d[which(d$aa_type=="N"),9], col="green", pch=16)
#points(d[which(d$aa_type=="F"),7], d[which(d$aa_type=="F"),9], col="orange", pch=16)
#abline(a=0, b=1, lty=3, lwd=3)
#title(main=infile, xlab="% of mutations, Unselected", ylab="% of mutations, Selected")
#legend("topleft",legend=c("Nonsynonymous", "Synonymous", "Nonsense", "Frameshift"), col=c("blue","red","green","orange"), pch=16, box.col="grey25")

#~Histogram~~~~~~~~~~~~~~~~~~~~#
#hist(d$count,40, main=infile, xlab="log2 Enrichment", col="grey")
#abline(v=mean(d$count), col="red", lty=2, lwd=3)

#dev.off()

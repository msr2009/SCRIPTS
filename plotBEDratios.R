d <- read.delim(infile, sep='\t', header=T)

#filter out low read-counts
#d <- d[-which(d[,6] <= 100),]

pdf(paste(strsplit(infile, split=".txt")[[1]][1], ".plots.pdf", sep=""), width=12, height=8)

#~~Position plot~~~~~~~~~~~~~~~~~~~~#
plot.new()
plot.window(xlim=c(0,max(d$aa_position)), ylim=c(min(d$ratio),max(d$ratio)))
axis(1, at=(c(seq(0,max(d$aa_position),by=50), max(d$aa_position))))
axis(2, cex.axis=2)
#plot synonymous changes
points(d$aa_position[which(d$aa_type=="S")], d$ratio[which(d$aa_type=="S")], col="red", pch=16, cex=1.5)
#plot nonsynonymous changes
points(d$aa_position[which(d$aa_type=="NS")], d$ratio[which(d$aa_type=="NS")], col="blue", pch=16, cex=1.5)
#plot nonsense changes
points(d$aa_position[which(d$aa_type=="N")], d$ratio[which(d$aa_type=="N")], col="green", pch=16, cex=1.5)
#plot frameshifts
points(d$aa_position[which(d$aa_type=="F")], d$ratio[which(d$aa_type=="F")], col="orange", pch=16, cex=1.5)

#title(main=infile, xlab="AA position", ylab="log2 Enrichment", cex.lab=1.25)

#~~Boxplots~~~~~~~~~~~~~~~~~~~~~~#
boxplot(d$ratio ~ d$aa_type, col=c("orange", "green", "blue", "red"), main=infile, ylab="log2 Enrichment", xlab="Mutation Type")

#~Scatters~~~~~~~~~~~~~~~~~~~~~~~#
plot.new()
plot.window(xlim=c(0,max(d[,7])), ylim=c(0,max(d[,9])))
axis(1); axis(2)
points(d[which(d$aa_type=="S"),7], d[which(d$aa_type=="S"),9], col="red", pch=16)
points(d[which(d$aa_type=="NS"),7], d[which(d$aa_type=="NS"),9], col="blue", pch=16)
points(d[which(d$aa_type=="N"),7], d[which(d$aa_type=="N"),9], col="green", pch=16)
points(d[which(d$aa_type=="F"),7], d[which(d$aa_type=="F"),9], col="orange", pch=16)
abline(a=0, b=1, lty=3, lwd=3)
title(main=infile, xlab="% of mutations, Unselected", ylab="% of mutations, Selected")
legend("topleft",legend=c("Nonsynonymous", "Synonymous", "Nonsense", "Frameshift"), col=c("blue","red","green","orange"), pch=16, box.col="grey25")

#~Histogram~~~~~~~~~~~~~~~~~~~~#
hist(d$ratio,40, main=infile, xlab="log2 Enrichment", col="grey")
abline(v=mean(d$ratio), col="red", lty=2, lwd=3)

dev.off()

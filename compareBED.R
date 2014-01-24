args <- commandArgs(TRUE)
input <- args[1]
sel <- args[2]
l <- args[3]

dati <- read.delim(input, sep='\t')
dats <- read.delim(sel, sep='\t')

dati.freq <- tabulate(dati$start[as.numeric(dati$start) < l])/length(dati$start)
dats.freq <- tabulate(dats$start[as.numeric(dats$start) < l])/length(dats$start)
enrichments <- log(dats.freq/dati.freq,2)

pdf('test.pdf', width=12, height=8)
plot(dati.freq, dats.freq)
abline(a=0,b=1,lty=2, col='red', lwd=3)
hist(enrichments, 50)
plot(seq(dats.freq), enrichments)
dev.off()


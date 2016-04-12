setwd("~/Desktop/")

#function for selection
#updates population frequencies for each generation
selection <- function(p, f, t, s){
  condition <- 1
  for( i in seq(t) ){
    p[,i+1] <- p[,i]*f[,condition]
    p[,i+1] <- p[,i+1]/sum(p[,i+1])
    if( i%%s == 0 ){
      condition <- 3-condition
      ###can add simulation to define new frequencies after backdilution
    }  
  }
  return(p)
}

#create initial population matrix, based on number strains
n <- 4
total_gen <- 100
sel_gen <- 5
pop <- matrix(data=NA, ncol=total_gen+1, nrow=n) 
pop[,1] <- 1/n

#create matrix of strain fitnesses
#for right now, let's calculate fitness based on constant fitness
#costs in each condition. In the future, I think we'll take this 
#matrix as input (with fitnesses for each strain).
cost1 <- .2
cost2 <- .2

# #we have these strains:
# #0) stable w/ ligand, stable w/o
# #1) stable w/ ligand, unstable w/o
# #2) unstable w/ ligand, stable w/o
# #3) unstable w/ ligand, unstable w/o
fitness <- matrix(c(1, 1, 1-cost1, 1-cost1, 1-cost2, 1, 1-cost2, 1))
dim(fitness) <- c(4,2)

# fitness <- function(strains){
#   s <- rnorm(n)
#   s <- (1-.7)/(max(s)-min(s))*(s-min(s))+.7
#   s.lig <- rnorm(n)
#   s.lig <- (1-.7)/(max(s.lig)-min(s.lig))*(s.lig-min(s.lig))+.7
#   return(cbind(s.lig, s))
# }
# f <- fitness(n)

colors <- c("yellow", "blue", "cyan", "magenta")
#postscript('tmp.eps', colormodel="cmyk")
sel <- selection(pop, f, total_gen, 5)
barplot(sel, col=colors, ylab="population frequency")
plot(sel, seq(1,total_gen), col=rainbow(n), ylab="frequency", xlab="generation", lty=1, lwd=.5)
plot(fitness, xlab="his+ligand", ylab="canavanine", pch=16, col="blue" )
#dev.off()

plot(0, 0, xlim=c(0,total_gen), ylim=c(0,1), type="n")
for(i in seq(n)){
  lines(seq(0,total_gen), sel[i,], col=rainbow(n)[i])
}

plot(f)

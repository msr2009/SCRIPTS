setwd("~/Desktop/")

#function for selection
#updates population frequencies for each generation
selection <- function(p, f, t, s, popsize, canR){
  condition <- 1
  for( i in seq(t) ){
    p[,i+1] <- p[,i]*f[,condition]
    p[,i+1] <- p[,i+1]/sum(p[,i+1])
    if( i%%s == 0 ){  ###switch condition every s generations
      condition <- 3-condition
      ###can add simulation to define new frequencies after backdilution
    }
    #add CanR cells to selection each generation (@2e-7/generation)
    if( canR ){
      floor(2e-7 * popsize)
    }
  }
  return(p)
}

#create initial population matrix, based on number strains
n <- 200
total_gen <- 100
sel_gen <- 5
pop <- matrix(data=NA, ncol=total_gen+1, nrow=n) 
pop[,1] <- 1/n
 
#create matrix of strain fitnesses
#for right now, let's calculate fitness based on constant fitness
#costs in each condition. In the future, I think we'll take this 
#matrix as input (with fitnesses for each strain).
cost1 <- .42
cost2 <- .42
costlow <- .25
# #we have these strains:
# #0) stable w/ ligand, stable w/o
# #1) stable w/ ligand, unstable w/o
# #2) unstable w/ ligand, stable w/o
# #3) unstable w/ ligand, unstable w/o
#f <- matrix(c(1, 1, 1-cost1, 1-cost1, 1-cost2, 1, 1-cost2, 1))
#dim(f) <- c(4,2)

fitness <- function(strains){
  s <- runif(n)
  s <- (1-costlow)/(max(s)-min(s))*(s-min(s))+costlow
  s.lig <- runif(n)
  s.lig <- (1-costlow)/(max(s.lig)-min(s.lig))*(s.lig-min(s.lig))+costlow
  return(cbind(s.lig, s))
}
f <- fitness(n)
colors <- rainbow(n)

#colors <- c("yellow", "blue", "cyan", "magenta")
#postscript('tmp.eps', colormodel="cmyk")
sel <- selection(pop, f, total_gen, 5)
par(mfrow=c(1,2))
b <- barplot(sel, col=colors, ylab="population frequency")
#plot(sel, seq(1,total_gen), col=rainbow(n), ylab="frequency", xlab="generation", lty=1, lwd=.5)
plot(f, xlab="his +ligand", ylab="canavanine -ligand", pch=16, cex=.5, col=rainbow(n), xlim=c(costlow,1), ylim=c(costlow,1) )
w <- which(sel[,21]>.05)
points(f[w,], lwd=2)
#dev.off()

library(qtl)
data(map10)
map10 <- map10[1:2]
set.seed(8789993)
simcross <- sim.cross(map10, n.ind=250, type="bc",
                      model=rbind(c(1, 50, 1.5), c(2, 50, 0)))
simcross$pheno[,1] <- simcross$pheno[,1] + rnorm(nind(simcross), 0, 2*simcross$qtlgeno[,2])
simcross <- calc.genoprob(simcross)
out <- scanonevar(simcross)
summary(out, format="allpeaks")

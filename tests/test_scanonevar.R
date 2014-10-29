library(qtl)
data(map10)
map10 <- map10[1:2]
set.seed(8789993)
simcross <- sim.cross(map10, n.ind=125, type="bc",
                      model=rbind(c(1, 50, 1.5), c(2, 50, 0)))
simcross$pheno[,1] <- simcross$pheno[,1] + rnorm(nind(simcross), 0, 2*simcross$qtlgeno[,2])
simcross <- calc.genoprob(simcross)
out <- scanonevar(simcross,
                  tol=0.01)
summary(out, format="allpeaks")

####

data(fake.bc)
fake.bc <- fake.bc[1:2,1:150] # only chr 1 and 2, and first 100 individuals
fake.bc <- calc.genoprob(fake.bc, step=5)
out <- scanonevar(fake.bc,
                  tol=0.01)
summary(out, format="allpeaks")
covar <- fake.bc$pheno[,c("sex", "age")]
out <- scanonevar(fake.bc, mean_covar=covar, var_covar=covar,
                  tol=0.01)
summary(out, format="allpeaks")

#########Simulate a vQTL on Chromosome 1########

chromo=1
qtl.position=14 # 50 cM
N=nind(fake.bc)
a1<-fake.bc$geno[[chromo]]$prob[,,1]
y <- fake.bc$pheno$pheno1
y <- y + rnorm(N,0,exp(a1[,qtl.position]))
out <- scanonevar(fake.bc, y, mean_covar=covar, var_covar=covar)
summary(out, format="allpeaks")

out <- scanonevar(fake.bc, y, mean_covar=covar,
                  tol=0.01)
summary(out, format="allpeaks")

out <- scanonevar(fake.bc, y, var_covar=covar,
                  tol=0.01)
summary(out, format="allpeaks")

out <- scanonevar(fake.bc, y,
                  tol=0.01)
summary(out, format="allpeaks")

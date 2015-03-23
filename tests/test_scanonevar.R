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

#########Simulate a vQTL on Chromosome 10 ########

data(fake.f2)

fake.f2 <- calc.genoprob(fake.f2, step = 2)

N = nind(fake.f2)

marker.vals <- fake.f2$geno[[10]]$data[,3]
marker.vals[is.na(marker.vals)] <- 1

fake.f2$pheno$phenotype <- fake.f2$pheno$phenotype + rnorm(N, 0, exp(marker.vals))

out <- scanonevar(fake.f2, pheno.col = 'phenotype')
plot(out, lodcolumn = 1:2)

out <- scanonevar(fake.f2,
									mean_covar = fake.f2$pheno$sex,
									var_covar = fake.f2$pheno$sex)
plot(out, lodcolumn = 1:2)

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



######### Simulate an additive mQTL on Chr5
	# a dominance mQTL on Chr7
	# an additive vQTL on Chr12
	# and a dominance vQTL on Chr14  ########

data(fake.f2)

fake.f2 <- calc.genoprob(fake.f2, step = 2)

N = nind(fake.f2)

marker1.vals <- fake.f2$geno[[5]]$data[,3]
marker1.vals[is.na(marker1.vals)] <- 1

marker2.vals <- fake.f2$geno[[7]]$data[,4]
marker2.vals[is.na(marker2.vals)] <- 1

marker3.vals <- fake.f2$geno[[12]]$data[,3]
marker3.vals[is.na(marker3.vals)] <- 1

marker4.vals <- fake.f2$geno[[14]]$data[,3]
marker4.vals[is.na(marker4.vals)] <- 1


fake.f2$pheno$phenotype1 <-
	rnorm(n = N, 25, 1) +
	marker1.vals

fake.f2$pheno$phenotype2 <-
	rnorm(n = N, 25, 1) +
	3*(marker2.vals == 1)

fake.f2$pheno$phenotype3 <-
	rnorm(n = N, 25, 1) +
	rnorm(n = N, mean = 0, sd = 3*marker3.vals)

fake.f2$pheno$phenotype4 <-
	rnorm(n = N, 25, 1) +
	rnorm(n = N, mean = 0, sd = 3*(marker4.vals != 1))

out1 <- scanonevar(fake.f2,
									pheno.col = 'phenotype1',
									dom = TRUE,
									mean_covar = fake.f2$pheno$sex,
									var_covar = fake.f2$pheno$sex)

plot(out1, bandcol = 'gray')

out2 <- scanonevar(fake.f2,
									pheno.col = 'phenotype2',
									dom = TRUE,
									mean_covar = fake.f2$pheno$sex,
									var_covar = fake.f2$pheno$sex)

plot(out2, bandcol = 'gray')

out3 <- scanonevar(fake.f2,
									pheno.col = 'phenotype3',
									dom = TRUE,
									mean_covar = fake.f2$pheno$sex,
									var_covar = fake.f2$pheno$sex)

plot(out3, bandcol = 'gray')

out4 <- scanonevar(fake.f2,
									pheno.col = 'phenotype4',
									dom = TRUE,
									mean_covar = fake.f2$pheno$sex,
									var_covar = fake.f2$pheno$sex)

plot(out4, bandcol = 'gray')

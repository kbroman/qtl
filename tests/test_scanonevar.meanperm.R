######### Simulate mQTL and vQTL modeled with scanonevar

library(qtl)
set.seed(27517)

data(fake.f2)

N = nind(fake.f2)

# make an additive mQTL on Chr5
marker1.vals <- fake.f2$geno[[5]]$data[,3] - 2
marker1.vals[is.na(marker1.vals)] <- 0

fake.f2$pheno$phenotype1 <- rnorm(n = N, 25 + 5*marker1.vals, 3)


simple.scan <- scanone(fake.f2, pheno.col = 'phenotype1')




var.scan.perms <- scanonevar.meanperm(fake.f2,
																			pheno.col = 'phenotype1',
																			num.mean.perms = 3,
																			dom = FALSE,
																			mean_covar = fake.f2$pheno$sex)




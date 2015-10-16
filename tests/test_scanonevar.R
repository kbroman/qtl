######### Simulate mQTL and vQTL modeled with scanonevar

library(qtl)
set.seed(27599)

data(fake.f2)

# # artificially expand the dataset by 5 fold
# for (chr in 1:nchr(fake.f2)) {
#
# 	this.data <- fake.f2$geno[[chr]]$data
# 	fake.f2$geno[[chr]]$data <- rbind(this.data, this.data, this.data, this.data, this.data)
#
# 	this.errors <- fake.f2$geno[[chr]]$errors
# 	fake.f2$geno[[chr]]$errors <- rbind(this.errors, this.errors, this.errors, this.errors, this.errors)
# }
# fake.f2$pheno <- rbind(fake.f2$pheno, fake.f2$pheno, fake.f2$pheno, fake.f2$pheno, fake.f2$pheno)

fake.f2 <- calc.genoprob(fake.f2, step = 2)

N = nind(fake.f2)
fake.f2$pheno$sex <- rbinom(n = N, size = 1, prob = 0.5)
fake.f2$pheno$age <- rnorm(n = N, mean = 10, sd = 1)
fake.f2$pheno$phen1 <- rnorm(n = N, mean = 20, sd = 2)


margin.plot(cross = fake.f2,
            focal.phenotype.name = 'phenotype',
            marginal.phen.names = list('sex', 'age'),
            marginal.marker.names = 'D17M88')


scan1 <- scanone(cross = fake.f2, chr = c(15:19, 'X'), pheno.col = 'phen1')
plot(scan1, bandcol = 'gray')

varscan1 <- scanonevar(cross = fake.f2,
                       mean.formula = formula('phen1 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
                       var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
                       chrs = c(15:19, 'X'))

plot(varscan1, scanone.for.comparison = scan1)





# an additive mQTL on Chr19
marker2.name <- colnames(fake.f2$geno$`19`$data)[2]
marker2.vals <- fake.f2$geno$`19`$data[,marker2.name] - 2
marker2.vals[is.na(marker2.vals)] <- 0

fake.f2$pheno$phen2 <- rnorm(n = N, 25 + marker2.vals, 3)

varscan2 <- scanonevar(cross = fake.f2,
                       mean.formula = formula('phen2 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
                       var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
                       chrs = c(15:19, 'X'))

summary(varscan2)
plot(varscan2)


predictive.plot(cross = fake.f2,
                mean.formula = formula('phen2 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
                var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
                marker.name = marker2.name,
                phen.name = 'sex')






# an additive vQTL on chr 15
marker3.name <- colnames(fake.f2$geno$`19`$data)[2]
marker3.vals <- get.genotypes.by.marker.name(cross = fake.f2, marker.name = marker3.name, as.matrix = FALSE) - 2
marker3.vals[is.na(marker3.vals)] <- 0

fake.f2$pheno$phen3 <- rnorm(n = N, 25, exp(0.5*marker3.vals))

varscan3 <- scanonevar(cross = fake.f2,
                        mean.formula = formula('phen3 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
                        var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
                        chrs = c(15:19, 'X'))

summary(varscan3)
plot(varscan3)

predictive.plot(cross = fake.f2,
                mean.formula = formula('phen3 ~ age + sex*mean.QTL.add + sex*mean.QTL.dom'),
                var.formula = formula('~sex + age + sex*var.QTL.add + sex*var.QTL.dom'),
                marker.name = marker3.name,
                phen.name = 'sex')


perms <- scanonevar.perm(cross = fake.f2,
                         mean.formula = formula('phen3 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
                         var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
                         n.perms = 50,
                         chrs = c(15:19, 'X'))


varscan3b <- convert.scanonevar.to.empirical.ps(scan = varscan3, 
                                                null.scan.maxes = perms)

# if the effects are really strong, the empricical p value will underflow R's float
# and it's log will be -Inf, so we can't plot....maybe should replace 0's with .Machine$double.eps?
# not likely to come up in real work, so I'll leave it as an error-throwing case for now
plot(varscan3b)


# todo X chromosome tests




# varscan1a.perms <- scanonevar.perm(cross = fake.f2,
# 																	 pheno.name = 'phenotype1',
# 																	 chrs = 18:20,
# 																	 num.perms = 15,
# 																	 mean.covar.names = c('sex', 'age'),
# 																	 var.covar.names = 'sex',
# 																	 quiet = FALSE)
#
# saveRDS(object = varscan1a.perms, file = 'test_varscan_perms.RDS')
# varscan1a.perms <- readRDS('test_varscan_perms.RDS')
#
# varscan1a.emp.ps <-
# 	lods.to.emp.ps.scanonevar(scan = varscan1a,
# 														null.scan.maxes = varscan1a.perms)
#
# plot(varscan1a.emp.ps)
#
# fitplot.scanonevar(cross = fake.f2,
# 									 varscan = varscan1a,
# 									 name.marker.to.plot = 'DXM64')


# better test for x chr stuff
B6.C58.cross <- read.cross(format = 'csv',
													 file = '~/Dropbox (ValdarLab)/vQTL_reanalysis/F2_B6_C58_all.csv',
													 genotypes = c('B', 'H', 'C'),
													 na.strings = c('O'),
													 convertX = FALSE)

# B6.C58.cross$geno$`17`$data[,1] <- rbinom(n = nind(B6.C58.cross),
# 																					size = 3,
# 																					p = 0.5)
# B6.C58.cross$geno$`17`$data[52,1] <- NA
# B6.C58.cross <- calc.genoprob(cross = B6.C58.cross, step = 2.0)

first.x.marker <- B6.C58.cross$geno$`X`$data[,1]
first.x.marker[is.na(first.x.marker)] <- 1
first.18.marker <- B6.C58.cross$geno$`18`$data[,1]
first.18.marker[is.na(first.18.marker)] <- 1
B6.C58.cross$pheno$bigenic <- rnorm(n = length(first.x.marker),
																		mean = first.x.marker + first.18.marker)

# B6.C58.cross$pheno$sex <- as.numeric(B6.C58.cross$pheno$sex)

vs <- scanonevar(cross = B6.C58.cross,
								 pheno.name = 'TOTREAR',
								 mean.covar.names = 'sex',
								 var.covar.names = 'sex')

plot(vs)

vs.perms <- scanonevar.perm(cross = B6.C58.cross,
														pheno.name = 'TOTREAR',
														mean.covar.names = 'sex',
														var.covar.names = 'sex',
														num.perms = 15,
														quiet = FALSE)

# so <- scanone(cross = B6.C58.cross,
# 							chr = c(16:19, 'X'),
# 							pheno.col = 'TOTDIST',
# 							addcovar = as.numeric(B6.C58.cross$pheno$sex))
#
# vs <- ConcatScanoneOnScanoneVar(varscan = vs, scanone = so)
#
# plot(vs, so)

vs.emp.ps <- convert.varscan.to.empirical.ps(scan = vs,
																						 null.scan.maxes = vs.perms)

plot(vs.emp.ps)



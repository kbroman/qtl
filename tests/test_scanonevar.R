######### Simulate mQTL and vQTL modeled with scanonevar
setwd('tests')

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
fake.f2$pheno$sex = rbinom(n = N, size = 1, prob = 0.5)
fake.f2$pheno$age = rnorm(n = N, mean = 10, sd = 1)


scan1a <- scanone(cross = fake.f2, chr = 1:5, pheno.col = 'phenotype')

# an additive mQTL on Chr5
marker1.vals <- fake.f2$geno[[20]]$data[,2] - 2
marker1.vals[is.na(marker1.vals)] <- 0

fake.f2$pheno$phenotype1 <- rnorm(n = N, 25 + 3*marker1.vals, 3)

varscan1a <- scanonevar(cross = fake.f2,
												pheno.name = 'phenotype1',
												chrs = 16:20,
												mean.covar.names = c('sex', 'age'),
												var.covar.names = 'sex',
												return.effects = TRUE)

saveRDS(object = varscan1a, file = 'test_varscan.RDS')
varscan1a <- readRDS('test_varscan.RDS')

summary(varscan1a)
plot(varscan1a)

varscan1a.perms <- scanonevar.perm(cross = fake.f2,
																	 pheno.name = 'phenotype1',
																	 chrs = 18:20,
																	 num.perms = 15,
																	 mean.covar.names = c('sex', 'age'),
																	 var.covar.names = 'sex',
																	 quiet = FALSE)

saveRDS(object = varscan1a.perms, file = 'test_varscan_perms.RDS')
varscan1a.perms <- readRDS('test_varscan_perms.RDS')

varscan1a.emp.ps <-
	lods.to.emp.ps.scanonevar(scan = varscan1a,
														null.scan.maxes = varscan1a.perms)

plot(varscan1a.emp.ps)

fitplot.scanonevar(cross = fake.f2,
									 varscan = varscan1a,
									 name.marker.to.plot = 'DXM64')


# better test for x chr stuff
B6.C58.cross <- read.cross(format = 'csv',
													 file = '~/Dropbox (ValdarLab)/vQTL_reanalysis/F2_B6_C58_all.csv',
													 genotypes = c('B', 'H', 'C'),
													 na.strings = c('O'),
													 convertX = FALSE)

B6.C58.cross <- calc.genoprob(cross = B6.C58.cross, step = 2.0)

first.x.marker <- B6.C58.cross$geno$`X`$data[,1]
first.x.marker[is.na(first.x.marker)] <- 1
B6.C58.cross$pheno$onx <- rnorm(n = length(first.x.marker),
																mean = first.x.marker,
																sd = 0.5)

vs <- scanonevar(cross = B6.C58.cross,
								 pheno.name = 'TOTDIST',
								 chrs = 16:20,
								 mean.covar.names = 'sex',
								 var.covar.names = 'sex')

so <- scanone(cross = B6.C58.cross,
							chr = c(16:19, 'X'),
							pheno.col = 'TOTDIST',
							addcovar = as.numeric(B6.C58.cross$pheno$sex))

vs <- ConcatScanoneOnScanoneVar(varscan = vs, scanone = so)

plot(vs, so)

vs.perms <- scanonevar.perm(cross = B6.C58.cross,
														pheno.name = 'TOTDIST',
														chrs = 16:20,
														mean.covar.names = 'sex',
														var.covar.names = 'sex',
														num.perms = 15,
														quiet = FALSE)

vs.emp.ps <- lods.to.emp.ps.scanonevar(scan = vs,
																			 null.scan.maxes = vs.perms)

plot(vs.emp.ps)

fitplot.scanonevar(cross = B6.C58.cross,
									 name.of.marker.to.plot = "17.45.986",
									 varscan = vs)

fitplot.scanonevar(cross = B6.C58.cross,
									 name.of.marker.to.plot = "17.45.986",
									 varscan = vs.emp.ps)

fitplot.scanonevar(cross = B6.C58.cross,
									 name.of.marker.to.plot = "X.106.858",
									 varscan = vs)

fitplot.scanonevar(cross = B6.C58.cross,
									 name.of.marker.to.plot = "X.106.858",
									 varscan = vs.emp.ps)


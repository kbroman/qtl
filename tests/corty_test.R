library(qtl)

# read data
{
B6.C58.cross <- read.cross(format = 'csv',
													 file = file.path('~',
													 								 'Dropbox (ValdarLab)',
													 								 'vQTL_reanalysis',
													 								 'F2_B6_C58_all.csv'),
													 genotypes = c('B', 'H', 'C'),
													 na.strings = c('O'),
													 convertX = FALSE)
B6.C58.cross <- calc.genoprob(cross = B6.C58.cross,
															step = 2)

male.B6.C58.cross <- read.cross(format = 'csv',
													 file = file.path('~',
													 								 'Dropbox (ValdarLab)',
													 								 'vQTL_reanalysis',
													 								 'F2_B6_C58_males.csv'),
													 genotypes = c('B', 'H', 'C'),
													 na.strings = c('O'),
													 convertX = FALSE)
male.B6.C58.cross <- calc.genoprob(cross = male.B6.C58.cross,
																	 step = 2)

female.B6.C58.cross <- read.cross(format = 'csv',
																file = file.path('~',
																								 'Dropbox (ValdarLab)',
																								 'vQTL_reanalysis',
																								 'F2_B6_C58_female.csv'),
																genotypes = c('B', 'H', 'C'),
																na.strings = c('O'),
																convertX = FALSE)
female.B6.C58.cross <- calc.genoprob(cross = female.B6.C58.cross,
																		 step = 2)

crosses <- list(all = B6.C58.cross,
								males = male.B6.C58.cross,
								females = female.B6.C58.cross)
}

# scanone on all mice and on each sex
for (cross.idx in 1:length(crosses)) {

	name <- names(crosses[cross.idx])
	so <- scanone(cross = crosses[[cross.idx]],
								pheno.col = 'TOTREAR')
	sop <- scanone(cross = crosses[[cross.idx]],
								 pheno.col = 'TOTREAR',
								 n.perm = 200)
	thresh <- quantile(sop, 0.95)
	plot(so,
			 ylim = c(0, max(so$lod, thresh)),
			 bandcol = 'lightgray',
			 main = paste('scanone', name))
	abline(h = thresh)
}



vs <- scanonevar(cross = crosses[[cross.idx]],
								 pheno.name = 'TOTREAR',
								 mean.covar.names = NULL,
								 var.covar.names = NULL)

saveRDS(object = vs, file = 'tests/test_vs.RDS')
vs <- readRDS(file = 'tests/test_vs.RDS')
plot(vs)

vs.perms <- scanonevar.perm(cross = B6.C58.cross,
														pheno.name = 'TOTREAR',
														mean.covar.names = 'sex',
														var.covar.names = 'sex',
														num.perms = 3,
														num.processes = 3,
														quiet = FALSE)
saveRDS(object = vs.perms, file = 'tests/test_vs_perms.RDS')
vs.perms <- readRDS(file = 'tests/test_vs_perms.RDS')

vs.emp.ps <- convert.varscan.to.empirical.ps(scan = vs,
																						 null.scan.maxes = vs.perms)

plot(vs.emp.ps)
plot(vs.emp.ps, chrs = 2)


# vs.cond <- scanonevar(cross = B6.C58.cross,
# 								 pheno.name = 'TOTREAR',
# 								 mean.covar.names = c('sex', 'X2.65.484'),
# 								 var.covar.names = c('sex', 'X2.65.484'))
#
# saveRDS(object = vs.cond, file = 'tests/test_vs_conditional.RDS')
vs.cond <- readRDS(file = 'tests/test_vs_conditional.RDS')
plot(vs.cond)


vs.cond.perms <- scanonevar.perm(cross = B6.C58.cross,
														pheno.name = 'TOTREAR',
														mean.covar.names = c('sex', 'X2.65.484'),
														var.covar.names = c('sex', 'X2.65.484'),
														num.perms = 15,
														quiet = FALSE)
saveRDS(object = vs.cond.perms, file = 'tests/test_vs_cond_perms.RDS')
vs.cond.perms <- readRDS(file = 'tests/test_vs_cond_perms.RDS')

vs.cond.emp.ps <- convert.varscan.to.empirical.ps(scan = vs.cond,
																									null.scan.maxes = vs.cond.perms)
saveRDS(object = vs.cond.emp.ps, file = 'tests/test_vs_cond_emp_ps.RDS')
vs.cond.emp.ps <- readRDS(file = 'tests/test_vs_cond_emp_ps.RDS')

plot(vs.cond.emp.ps)


st <- scantwo(cross = B6.C58.cross, pheno.col = 'TOTREAR')


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

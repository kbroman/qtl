library(qtl)

# # read data
# {
# B6.C58.cross <- read.cross(format = 'csv',
# 													 file = file.path('rds',
# 													 								 'F2_B6_C58_all.csv'),
# 													 genotypes = c('B', 'H', 'C'),
# 													 na.strings = c('O'),
# 													 convertX = FALSE)
# B6.C58.cross <- calc.genoprob(cross = B6.C58.cross,
# 															step = 2)
# 
# male.B6.C58.cross <- read.cross(format = 'csv',
# 													 file = file.path('~',
# 													 								 'Dropbox (ValdarLab)',
# 													 								 'vQTL_reanalysis',
# 													 								 'F2_B6_C58_males.csv'),
# 													 genotypes = c('B', 'H', 'C'),
# 													 na.strings = c('O'),
# 													 convertX = FALSE)
# male.B6.C58.cross <- calc.genoprob(cross = male.B6.C58.cross,
# 																	 step = 2)
# 
# female.B6.C58.cross <- read.cross(format = 'csv',
# 																file = file.path('~',
# 																								 'Dropbox (ValdarLab)',
# 																								 'vQTL_reanalysis',
# 																								 'F2_B6_C58_female.csv'),
# 																genotypes = c('B', 'H', 'C'),
# 																na.strings = c('O'),
# 																convertX = FALSE)
# female.B6.C58.cross <- calc.genoprob(cross = female.B6.C58.cross,
# 																		 step = 2)
# 
# crosses <- list(all = B6.C58.cross,
# 								males = male.B6.C58.cross,
# 								females = female.B6.C58.cross)
# }

B6.C58.cross <- readRDS('rds/B6_C58_cross.RDS')

margin.plot(cross = B6.C58.cross, 
            focal.phenotype.name = 'TOTDIST',
            marginal.phen.names = list('sex', 'AMBEPIS'), 
            marginal.marker.names = list("X1.52.903"))

B6.C58.cross$pheno$sqrt.totrear <- sqrt(B6.C58.cross$pheno$TOTREAR)

so <- scanone(cross = B6.C58.cross,
              chr = 1:3,
              pheno.col = 'sqrt.totrear')

a <- newScanonevar(cross = B6.C58.cross,
                   mean.formula = sqrt.totrear ~ sex + mean.QTL.add + mean.QTL.dom,
                   var.formula = ~sex + var.QTL.add + var.QTL.dom,
                   return.effects = TRUE,
                   return.effect.ses = TRUE,
                   return.effect.ps = TRUE,
                   chrs = c(1, 2, 3), 
                   exclusion.window = 0.8)

plot(a, scanone.for.comparison = so)

predictive.plot(varscan = a, cross = B6.C58.cross,
                phen.name = 'sex', marker.name = 'X1.58.218')

summary(a, thresh = 3)

fitplot.scanonevar(cross = B6.C58.cross,
                   name.of.marker.to.plot = 'chr01_X1.63.544',
                   varscan = a,
                   col = B6.C58.cross$pheno$sex,
                   sample.based = TRUE,
                   phenotype.name = 'sqrt.totrear')

levels(B6.C58.cross$pheno$sex) <- c('red', 'blue')
fitplot.scanonevar(cross = B6.C58.cross,
                   name.of.marker.to.plot = 'chr01_X1.63.544',
                   varscan = a,
                   col = B6.C58.cross$pheno$sex,
                   fit.based = TRUE)




# working on wald test
theta <- both.alt.fit$coef[3:4]
R <- diag(2)
r <- c(0, 0)
V <- summary(both.alt.fit, correlation = TRUE)$cov.scaled[3:4, 3:4]
t(R %*% theta - r) %*% solve(R %*% (V) %*% t(R)) %*% (R %*% theta - r)

fitplot.scanonevar(cross = B6.C58.cross, 
                   name.of.marker.to.plot = 'chr01_X1.58.218',
                   varscan = a)

b <- newScanonevar.perm(cross = B6.C58.cross,
                        mean.formula = sqrt.totrear ~ sex + mean.QTL.add + mean.QTL.dom,
                        var.formula = ~sex + var.QTL.add + var.QTL.dom + X3.58.232,
                        chrs = 1,
                        n.perms = 10)
b

c <- newConvenientScanonevar(cross = B6.C58.cross,
                             mean.formula = sqrt.totrear ~ sex*(QTL.add + QTL.dom) + X1.152.072,
                             var.formula = ~sex ,
                             n.perms = 24,
                             n.cores = 8)
plot(c)
 
# sov <- scanonevar(cross = B6.C58.cross,
#                   pheno.name = 'sqrt.totrear',
# #                   chrs = 18:20,
#                   mean.covar.names = c('sex', 'X2.136.176', 'X7.32.784'),
#                   var.covar.names = 'sex')

plot(sov)


svp <- scanonevar.perm(cross = B6.C58.cross,
                       pheno.name = 'sqrt.totrear',
#                        chrs = 15:20,
                       mean.covar.names = c('sex', 'X2.136.176', 'X7.32.784'),
                       var.covar.names = 'sex',
                       num.perms = 4)


vqtl_ConvertLODsToEmpPs(scan = sov,
                        null.scan.maxes = svp)

vseps <- ConvenientScanoneVar(cross = B6.C58.cross,
                              pheno.name = 'sqrt.totrear',
                              #chrs = 15:20,
                              mean.covar.names = c('sex', 'X2.136.176', 'X7.32.784'),
                              var.covar.names = 'sex',
                              num.perms = 50,
                              quiet = FALSE)

plot(vseps)

summary(vseps)

fitplot.scanonevar(cross = B6.C58.cross,
                   name.of.marker.to.plot = 'X16.56.11',
                   varscan = vseps)

plot(vseps, chr = 2)
plot(vseps, chr = 7)

# fitplots for full model hits
fitplot.scanonevar(cross = B6.C58.cross,
                   name.of.marker.to.plot = 'X2.65.484',
                   varscan = vseps)
fitplot.scanonevar(cross = B6.C58.cross,
                   name.of.marker.to.plot = 'X2.120.543',
                   varscan = vseps)
fitplot.scanonevar(cross = B6.C58.cross,
                   name.of.marker.to.plot = 'X7.32.784',
                   varscan = vseps)

# fitplots for mean model hit
fitplot.scanonevar(cross = B6.C58.cross,
                   name.of.marker.to.plot = 'X2.136.176',
                   varscan = vseps)
fitplot.scanonevar(cross = B6.C58.cross,
                   name.of.marker.to.plot = 'X7.32.784',
                   varscan = vseps)

# fitplots for var model hits
fitplot.scanonevar(cross = B6.C58.cross,
                   name.of.marker.to.plot = 'X2.65.484',
                   varscan = vseps)



cond.vseps <- ConvenientScanoneVar(cross = B6.C58.cross,
                                   pheno.name = 'sqrt.totrear',
                                   mean.covar.names = c('sex', 'X2.65.484'),
                                   num.perms = 52,
                                   num.processes = 4)

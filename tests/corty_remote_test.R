library(qtl)
library(CortyKit)

B6.C58.cross <- EnsureExists(var = B6.C58.cross,
                             expr = {B6.C58.cross <- read.cross(format = 'csv',
                                                file = file.path('~',
                                                                 'Dropbox (ValdarLab)',
                                                                 'vQTL_reanalysis',
                                                                 'F2_B6_C58_all.csv'),
                                                genotypes = c('B', 'H', 'C'),
                                                na.strings = c('O'),
                                                convertX = FALSE)
                               B6.C58.cross <- calc.genoprob(cross = B6.C58.cross,
                                                             step = 2) },
                             dir = 'rds')

B6.C58.cross$pheno$sqrt.totrear <- sqrt(B6.C58.cross$pheno$TOTREAR)
B6.C58.cross$pheno$rint.totrear <- rint(B6.C58.cross$pheno$TOTREAR)


sov6 <- EnsureExists(var = sov6,
                     expr = newConvenientScanonevar(cross = B6.C58.cross,
                                                    mean.formula = sqrt(TOTREAR) ~ sex*(QTL.add + QTL.dom) + X2.65.484,
                                                    var.formula = ~sex + QTL.add + QTL.dom + X2.65.484,
                                                    n.perms = 8*10,
                                                    n.cores = 8),
                     dir = 'rds')





sov4 <- EnsureExists(var = sov4,
                     expr = newConvenientScanonevar(cross = B6.C58.cross,
                                                    mean.formula = sqrt(TOTREAR) ~ sex*(QTL.add + QTL.dom) + X2.65.484,
                                                    var.formula = ~sex + QTL.add + QTL.dom,
                                                    n.perms = 8*15,
                                                    n.cores = 8),
                     dir = 'rds')
plot(sov4)


sov <- scanonevar(cross = B6.C58.cross,
                  pheno.name = 'sqrt.totrear',
                  chrs = 15:20,
                  var.add = FALSE,
                  var.dom = FALSE,
                  mean.covar.names = 'sex',
                  var.covar.names = 'sex')
                                 

vsperms <- scanonevar.perm(cross = B6.C58.cross, 
                           pheno.name = 'sqrt.totrear', 
                           var.add = FALSE, 
                           var.dom = FALSE, 
                           mean.covar.names = 'sex', 
                           var.covar.names = 'sex', 
                           num.perms = 16*4,
                           num.processes = 16)


eps <- vqtl_ConvertLODsToEmpPs(scan = sov, null.scan.maxes = vsperms)


vseps <- ConvenientScanoneVar(cross = B6.C58.cross,
                              pheno.name = 'sqrt.totrear',
#                               chrs = 15:20,
                              var.add = FALSE,
                              var.dom = FALSE,
                              mean.covar.names = 'sex',
                              var.covar.names = 'sex',
                              num.perms = 16*6,
                              num.processes = 16)

plot(vseps)

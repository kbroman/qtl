######### Simulate mQTL and vQTL modeled with scanonevar

#library(qtl)
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




# an additive mQTL on Chr5

marker1.vals <- fake.f2$geno[[5]]$data[,3] - 2
marker1.vals[is.na(marker1.vals)] <- 0

fake.f2$pheno$phenotype1 <- rnorm(n = N, 25 + 5*marker1.vals, 3)

var.scan1a <- scanonevar(fake.f2,
												pheno.col = 'phenotype1',
												chrs = 4:6,
												dom = FALSE)

plot(var.scan1a, bandcol = 'gray')
fitplot.scanonevar(cross = fake.f2, var.scan = var.scan1a, marker.name = 'D5M391')

# df <- data.frame(phen1 = fake.f2$pheno$phenotype1, marker1 = marker1.vals)
#
# dg1 <- dglm(formula = phen1 ~ marker1,
# 						dformula = ~ 1,
# 						data = df,
# 						family = gaussian)
# ln.lik <- -0.5*dg1$m2loglik
# log10.lik <- ln.lik / log(10)

# var.scan1b <- scanonevar(fake.f2,
# 												 pheno.col = 'phenotype1',
# # 												 use.dglm.package = FALSE,
# # 												 use.custom.em = TRUE,
# 												 chrs = 4:6,
# 												 dom = TRUE)
#
# plot(var.scan1b, bandcol = 'gray')
# fitplot.scanonevar(cross = fake.f2, var.scan = var.scan1b, marker.name = 'D5M391')





# a dominance mQTL on Chr7
marker2.vals <- fake.f2$geno[[7]]$data[,4]-2
marker2.vals[is.na(marker2.vals)] <- 0
marker2.dom <- marker2.vals == 0

fake.f2$pheno$phenotype2 <- rnorm(n = N, 25 + 5*(marker2.vals != -1), 3)

var.scan2a <- scanonevar(fake.f2,
												 pheno.col = 'phenotype2',
												 dom = TRUE,
# 												 use.dglm.package = TRUE,
# 												 use.custom.em = FALSE,
												 chrs = 6:8)

plot(var.scan2a, bandcol = 'gray')
fitplot.scanonevar(cross = fake.f2, var.scan = var.scan2a, marker.name = 'D7M285')

# var.scan2b <- scanonevar(fake.f2,
# 												 pheno.col = 'phenotype2',
# 												 dom = TRUE,
# # 												 use.dglm.package = FALSE,
# # 												 use.custom.em = TRUE,
# 												 chrs = 6:8)
#
# #plot(var.scan2B, bandcol = 'gray')
# fitplot.scanonevar(cross = fake.f2, var.scan = var.scan2b, marker.name = 'D7M285')






# an additive vQTL on Chr12
marker3.vals <- fake.f2$geno[[12]]$data[,3] - 2
marker3.vals[is.na(marker3.vals)] <- - 1
marker3.dom <- marker3.vals == 0

fake.f2$pheno$phenotype3 <- rnorm(n = N, 25, sd = 3*(marker3.vals + 2))

var.scan3a <- scanonevar(fake.f2,
												 pheno.col = 'phenotype3',
												 dom = FALSE,
# 												 use.dglm.package = TRUE,
# 												 use.custom.em = FALSE,
												 chrs = 11:13)

plot(var.scan3a, bandcol = 'gray', legend.pos = c('topleft', 'bottomleft'))
fitplot.scanonevar(cross = fake.f2, var.scan = var.scan3a, marker.name = 'D12M52')

# var.scan3b <- scanonevar(fake.f2,
# 												 pheno.col = 'phenotype3',
# 												 dom = FALSE,
# # 												 use.dglm.package = FALSE,
# # 												 use.custom.em = TRUE,
# 												 chrs = 11:13)
#
# #plot(var.scan3b, bandcol = 'gray', legend.pos = c('topleft', 'bottomleft'))
# fitplot.scanonevar(cross = fake.f2, var.scan = var.scan3b, marker.name = 'D12M52')




# a dominance vQTL on Chr14  ########

marker4.vals <- fake.f2$geno[[14]]$data[,3] - 2
marker4.vals[is.na(marker4.vals)] <- -1
marker4.dom <- marker4.vals == 0

fake.f2$pheno$phenotype4 <-	rnorm(n = N,
																	25,
																	sd = 3*((marker4.vals != -1)+1))

var.scan4a <- scanonevar(fake.f2,
												 pheno.col = 'phenotype4',
# 												 use.dglm.package = TRUE,
# 												 use.custom.em = FALSE,
												 dom = TRUE,
												 chrs = 13:15)

plot(var.scan4a, bandcol = 'gray', legend.pos = 'topleft')
fitplot.scanonevar(cross = fake.f2, var.scan = var.scan4a, marker.name = 'D14M160')

# var.scan4b <- scanonevar(fake.f2,
# 												 pheno.col = 'phenotype4',
# 												 use.dglm.package = FALSE,
# 												 use.custom.em = TRUE,
# 												 dom = TRUE,
# 												 chrs = 13:15)
#
# #plot(var.scan4b, bandcol = 'gray', legend.pos = 'topleft')
# fitplot.scanonevar(cross = fake.f2, var.scan = var.scan4b, marker.name = 'D14M160')



# a dominance mQTL and dominance vQTL on Chr15  ########

marker5.vals <- fake.f2$geno[[15]]$data[,1] - 2
marker5.vals[is.na(marker5.vals)] <- 0
marker5.dom <- marker5.vals == 0

fake.f2$pheno$phenotype5 <-	rnorm(n = N,
																	mean = 25 + 5*(marker5.vals != -1),
																	sd = 3*((marker5.vals != -1)+1))

var.scan5a <- scanonevar(fake.f2,
												 pheno.col = 'phenotype5',
# 												 use.dglm.package = TRUE,
# 												 use.custom.em = FALSE,
												 dom = TRUE,
												 chrs = 14:16)

plot(var.scan5a, bandcol = 'gray', legend.pos = 'topleft')
fitplot.scanonevar(cross = fake.f2, var.scan = var.scan5a, marker.name = 'D15M5')

# var.scan5b <- scanonevar(fake.f2,
# 												 pheno.col = 'phenotype5',
# # 												 use.dglm.package = FALSE,
# # 												 use.custom.em = TRUE,
# 												 dom = TRUE,
# 												 chrs = 14:16)
#
# #plot(var.scan5b, bandcol = 'gray', legend.pos = 'topleft')
# fitplot.scanonevar(cross = fake.f2, var.scan = var.scan5b, marker.name = 'D15M5')


# a dominance vQTL on Chr15 (same marker as before, but no mean QTL) ########

fake.f2$pheno$phenotype6 <-	rnorm(n = N,
																	mean = 25,
																	sd = 3*((marker5.vals != -1)+1))

var.scan6 <- scanonevar(fake.f2,
												pheno.col = 'phenotype6',
												dom = TRUE,
												chrs = 14:16)

plot(var.scan6, bandcol = 'gray', legend.pos = 'topleft')
fitplot.scanonevar(cross = fake.f2, var.scan = var.scan6, marker.name = 'D15M5')

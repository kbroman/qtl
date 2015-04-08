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
plot(simple.scan)



var.scan.perms <- scanonevar.varperm(fake.f2,
																		 pheno.col = 'phenotype1',
																		 num.var.perms = 2,
																		 dom = FALSE,
																		 mean_covar = fake.f2$pheno$sex)


library(dplyr)

meanperms <- rbind_all(var.scan.perms)



meanperms %>% group_by(perm.num) %>% summarise(full = max(lod.full), mean = max(lod.mean))


###### Timing  ##########

scan1 <- system.time(var.scan <- scanonevar(fake.f2,
																						pheno.col = 'phenotype1',
																						dom = FALSE))[3]


scan5 <- system.time(var.scan.perms5 <- scanonevar.meanperm(fake.f2,
																														pheno.col = 'phenotype1',
																														num.mean.perms = 3,
																														dom = FALSE))[3]


scan10 <- system.time(var.scan.perms10 <- scanonevar.meanperm(fake.f2,
																															pheno.col = 'phenotype1',
																															num.mean.perms = 10,
																															dom = FALSE))[3]


scan15 <- system.time(var.scan.perms15 <- scanonevar.meanperm(fake.f2,
																															pheno.col = 'phenotype1',
																															num.mean.perms = 15,
																															dom = FALSE))[3]


scan20 <- system.time(var.scan.perms20 <- scanonevar.meanperm(fake.f2,
																															pheno.col = 'phenotype1',
																															num.mean.perms = 20,
																															dom = FALSE))[3]

# looks like about 45 s per permutation
# so 1000 perms will take...750 mins = 13 hours compute time
# divide that over 100 computers.....
# ten perms per computer....
# 7 minutes per computer...



perm.maxes <- meanperms %>% group_by(perm.num) %>% summarise(max.lod.full = max(lod.full),
																														 max.lod.mean = max(lod.mean),
																														 max.lod.disp = max(lod.disp))


hist(perm.maxes$max.lod.full)
hist(perm.maxes$max.lod.mean, add = TRUE)
hist(perm.maxes$max.lod.disp, add = TRUE)

plot(var.scan)

library(stringr)

meanperms$chr <- str_pad(string = meanperms$chr, width = 2, side = 'left', pad = '0')
test <- meanperms %>% group_by(chr, pos) %>% summarise(lod = quantile(lod.full, probs = 1))
class(test) <- c('scanone', 'data.frame')
plot(test, bandcol = 'gray')


#### timing with dominance effects ####

scan1 <- system.time(var.scan <- scanonevar(fake.f2,
																						pheno.col = 'phenotype1',
																						dom = TRUE))[3]


scan5 <- system.time(var.scan.perms5 <- scanonevar.meanperm(fake.f2,
																														pheno.col = 'phenotype1',
																														num.mean.perms = 5,
																														dom = TRUE))[3]


scan10 <- system.time(var.scan.perms10 <- scanonevar.meanperm(fake.f2,
																															pheno.col = 'phenotype1',
																															num.mean.perms = 10,
																															dom = TRUE))[3]


##############################################################
# R code for "A brief tour of R/qtl"
#
# Karl W Broman, kbroman.wisc.edu
# University of Wisconsin Madison
#
# http://www.rqtl.org
#
# 21 March 2012
##############################################################

save.image()

library(qtl)

ls()

help(read.cross)
?read.cross

############################################################
# Example 1: Hypertension
############################################################
data(hyper)
ls()
?hyper

summary(hyper)

nind(hyper)
nphe(hyper)
nchr(hyper)
totmar(hyper)
nmar(hyper)

plot(hyper)

plotMissing(hyper)
plotMap(hyper)
plotPheno(hyper, pheno.col=1)

plotMap(hyper, chr=c(1, 4, 6, 7, 15), show.marker.names=TRUE)

plotMissing(hyper, reorder=TRUE)

hyper <- drop.nullmarkers(hyper)
totmar(hyper)

hyper <- est.rf(hyper)
plotRF(hyper)
plotRF(hyper, chr=c(1,4))

plotRF(hyper, chr=6)
plotMissing(hyper, chr=6)

newmap <- est.map(hyper, error.prob=0.01)
plotMap(hyper, newmap)

hyper <- replace.map(hyper, newmap)

hyper <- calc.errorlod(hyper, error.prob=0.01)

top.errorlod(hyper)

plotGeno(hyper, chr=16, ind=c(24:34, 71:81))

plotInfo(hyper)
plotInfo(hyper, chr=c(1,4,15))
plotInfo(hyper, chr=c(1,4,15), method="entropy")
plotInfo(hyper, chr=c(1,4,15), method="variance")

hyper <- calc.genoprob(hyper, step=1, error.prob=0.01)

out.em <- scanone(hyper)
out.hk <- scanone(hyper, method="hk")

hyper <- sim.geno(hyper, step=2, n.draws=16, error.prob=0.01)
out.imp <- scanone(hyper, method="imp")

summary(out.em)
summary(out.em, threshold=3)
summary(out.hk, threshold=3)
summary(out.imp, threshold=3)

max(out.em)
max(out.hk)
max(out.imp)

plot(out.em, chr=c(1,4,15))
plot(out.em, out.hk, out.imp, chr=c(1,4,15))
plot(out.em, chr=c(1,4,15))
plot(out.hk, chr=c(1,4,15), col="blue", add=TRUE)
plot(out.imp, chr=c(1,4,15), col="red", add=TRUE)

operm.hk <- scanone(hyper, method="hk", n.perm=1000)

summary(operm.hk, alpha=0.05)

summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE)

save.image()

hyper <- calc.genoprob(hyper, step=5, error.prob=0.01)

out2.hk <- scantwo(hyper, method="hk")

summary(out2.hk, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6))

summary(out2.hk, thresholds=c(6.0, 4.7, Inf, 4.7, 2.6))

plot(out2.hk)
plot(out2.hk, chr=c(1,4,6,15))

max(out2.hk)

operm2.hk <- scantwo(hyper, method="hk", n.perm=100)

summary(operm2.hk)

summary(out2.hk, perms=operm2.hk, pvalues=TRUE,
        alphas=c(0.05, 0.05, 0, 0.05, 0.05))

chr <- c(1, 1, 4, 6, 15)
pos <- c(50, 76, 30, 70, 20)
qtl <- makeqtl(hyper, chr, pos)

my.formula <- y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q4:Q5
out.fitqtl <- fitqtl(hyper, qtl=qtl, formula=my.formula)
summary(out.fitqtl)

ls()
rm(list=ls())

############################################################
# Example 2: Genetic mapping
############################################################
data(badorder)
summary(badorder)
plot(badorder)

badorder <- est.rf(badorder)
plotRF(badorder)

plotRF(badorder, chr=1)

newmap <- est.map(badorder, verbose=TRUE)
plotMap(badorder, newmap)

plotRF(badorder, chr=2:3)

pull.map(badorder, chr=2)
pull.map(badorder, chr=3)

badorder <- movemarker(badorder, "D2M937", 3, 48)
badorder <- movemarker(badorder, "D3M160", 2, 28.8)

plotRF(badorder, chr=2:3)

rip1 <- ripple(badorder, chr=1, window=6)
summary(rip1)

rip2 <- ripple(badorder, chr=1, window=3, err=0.01, method="likelihood")
summary(rip2)

badorder.rev <- switch.order(badorder, 1, rip1[2,])
rip1r <- ripple(badorder.rev, chr=1, window=6)
summary(rip1r)

badorder.rev <- switch.order(badorder.rev, 1, rip1r[2,])
rip2r <- ripple(badorder.rev, chr=1, window=3, err=0.01)
summary(rip2r)

badorder.rev <- est.rf(badorder.rev)
plotRF(badorder.rev, 1)

############################################################
# Example 3: Listeria susceptibility
############################################################
data(listeria)
summary(listeria)
plot(listeria)
plotMissing(listeria)

listeria$pheno$logSurv <- log(listeria$pheno[,1])
plot(listeria)

listeria <- est.rf(listeria)
plotRF(listeria)
plotRF(listeria, chr=c(5,13))

newmap <- est.map(listeria, error.prob=0.01)
plotMap(listeria, newmap)
listeria <- replace.map(listeria, newmap)

listeria <- calc.errorlod(listeria, error.prob=0.01)
top.errorlod(listeria)
top.errorlod(listeria, cutoff=3.5)
plotGeno(listeria, chr=13, ind=61:70, cutoff=3.5)

listeria <- calc.genoprob(listeria, step=2)
out.2p <- scanone(listeria, pheno.col=3, model="2part", upper=TRUE)

summary(out.2p)
summary(out.2p, threshold=4.5)

summary(out.2p, format="allpeaks", threshold=3)
summary(out.2p, format="allpeaks", threshold=c(4.5,3,3))

plot(out.2p)
plot(out.2p, lodcolumn=2)
plot(out.2p, lodcolumn=1:3, chr=c(1,5,13,15))

operm.2p <- scanone(listeria, model="2part", pheno.col=3,
                    upper=TRUE, n.perm=25)
summary(operm.2p, alpha=0.05)

summary(out.2p, format="allpeaks", perms=operm.2p,
        alpha=0.05, pvalues=TRUE)

y <- listeria$pheno$logSurv
my <- max(y, na.rm=TRUE)
z <- as.numeric(y==my)
y[y==my] <- NA
listeria$pheno$logSurv2 <- y
listeria$pheno$binary <- z
plot(listeria)

out.mu <- scanone(listeria, pheno.col=4)
plot(out.mu, out.2p, lodcolumn=c(1,3), chr=c(1,5,13,15), col=c("blue","red"))

out.p <- scanone(listeria, pheno.col=5, model="binary")
plot(out.p, out.2p, lodcolumn=c(1,2), chr=c(1,5,13,15), col=c("blue","red"))

out.p.alt <- scanone(listeria, pheno.col=as.numeric(listeria$pheno$T264==264),
                     model="binary")

out.np1 <- scanone(listeria, model="np", ties.random=TRUE)
out.np2 <- scanone(listeria, model="np", ties.random=FALSE)

plot(out.np1, out.np2, col=c("blue","red"))
plot(out.2p, out.np1, out.np2, chr=c(1,5,13,15))

############################################################
# Example 4: Covariates in QTL mapping
############################################################
data(fake.bc)
summary(fake.bc)
plot(fake.bc)

fake.bc <- calc.genoprob(fake.bc, step=2.5)
out.nocovar <- scanone(fake.bc, pheno.col=1:2)

sex <- fake.bc$pheno$sex
out.acovar <- scanone(fake.bc, pheno.col=1:2, addcovar=sex)

summary(out.nocovar, threshold=3, format="allpeaks")
summary(out.acovar, threshold=3, format="allpeaks")

plot(out.nocovar, out.acovar, chr=c(2, 5))
plot(out.nocovar, out.acovar, chr=c(2, 5), lodcolumn=2)

out.icovar <- scanone(fake.bc, pheno.col=1:2, addcovar=sex, intcovar=sex)

summary(out.icovar, threshold=3, format="allpeaks")

plot(out.acovar, out.icovar, chr=c(2,5), col=c("blue", "red"))
plot(out.acovar, out.icovar, chr=c(2,5), lodcolumn=2,
     col=c("blue", "red"))

out.sexint <- out.icovar - out.acovar
plot(out.sexint, lodcolumn=1:2, chr=c(2,5), col=c("green", "purple"))

seed <- ceiling(runif(1, 0, 10^8))
set.seed(seed)
operm.acovar <- scanone(fake.bc, pheno.col=1:2, addcovar=sex,
                        method="hk", n.perm=100)
set.seed(seed)
operm.icovar <- scanone(fake.bc, pheno.col=1:2, addcovar=sex,
                        intcovar=sex, method="hk", n.perm=100)

operm.sexint <- operm.icovar - operm.acovar

summary(operm.sexint, alpha=c(0.05, 0.20))

summary(out.sexint, perms=operm.sexint, alpha=0.1,
        format="allpeaks", pvalues=TRUE)

############################################################
# Example 5: Multiple QTL mapping
############################################################
rm(list=ls())
data(hyper)

hyper <- sim.geno(hyper, step=2.5, n.draws=16, err=0.01)

out1 <- scanone(hyper, method="imp")
plot(out1)

max(out1)

find.marker(hyper, 4, 29.5)

g <- pull.geno(hyper)[,"D4Mit164"]
mean(is.na(g))

g <- pull.geno(fill.geno(hyper))[,"D4Mit164"]

out1.c4 <- scanone(hyper, method="imp", addcovar=g)

plot(out1, out1.c4, col=c("blue", "red"))

plot(out1.c4 - out1, ylim=c(-3,3))
abline(h=0, lty=2, col="gray")

out1.c4i <- scanone(hyper, method="imp", addcovar=g, intcovar=g)

plot(out1.c4i - out1.c4)

out2 <- scantwo(hyper, method="imp")

summary(out2, thr=c(6.0, 4.7, Inf, 4.7, 2.6))

summary( subset(out2, chr=1) )

summary( subset(out2, chr=c(7,15)) )

plot(out2, chr=c(1,4,6,7,15))

plot(out2, chr=1, lower="cond-add")
plot(out2, chr=c(6,15), lower="cond-int")
plot(out2, chr=c(7,15), lower="cond-int")

out2.c4 <- scantwo(hyper, method="imp", addcovar=g, chr=c(1,6,7,15))

summary(out2.c4, thr=c(6.0, 4.7, Inf, 4.7, 2.6))
summary( subset(out2.c4, chr=1) )
summary( subset(out2.c4, chr=c(7,15)) )

plot(out2.c4)
plot(out2.c4, chr=1, lower="cond-int")
plot(out2.c4, chr=c(6,15), lower="cond-int")
plot(out2.c4, chr=c(7,15), lower="cond-int")

out2sub <- subset(out2, chr=c(1,6,7,15))
plot(out2.c4 - out2sub, allow.neg=TRUE, lower="cond-int")

qc <- c(1, 1, 4, 6, 15)
qp <- c(43.3, 78.3, 30.0, 62.5, 18.0)
qtl <- makeqtl(hyper, chr=qc, pos=qp)

myformula <- y ~ Q1+Q2+Q3+Q4+Q5 + Q4:Q5

out.fq <- fitqtl(hyper, qtl=qtl, formula = myformula)
summary(out.fq)

out.fq <- fitqtl(hyper, qtl=qtl, formula = myformula, drop=FALSE, get.ests=TRUE)
summary(out.fq)

revqtl <- refineqtl(hyper, qtl=qtl, formula = myformula)

revqtl

plot(revqtl)

out.fq2 <- fitqtl(hyper, qtl=revqtl, formula=myformula)
summary(out.fq2)

out1.c4r <- addqtl(hyper, qtl=revqtl, formula=y~Q3)

plot(out1.c4, out1.c4r, col=c("blue", "red"))

plot(out1.c4r - out1.c4, ylim=c(-1.7, 1.7))
abline(h=0, lty=2, col="gray")

out2.c4r <- addpair(hyper, qtl=revqtl, formula=y~Q3, chr=c(1,6,7,15))

plot(out2.c4r - out2.c4, lower="cond-int", allow.neg=TRUE)

out.1more <- addqtl(hyper, qtl=revqtl, formula=myformula)
plot(out.1more)

out.iw4 <- addqtl(hyper, qtl=revqtl, formula=y~Q1+Q2+Q3+Q4+Q5+Q4:Q5+Q6+Q5:Q6)
plot(out.iw4)

out.2more <- addpair(hyper, qtl=revqtl, formula=myformula, chr=c(2,5,7,15))

plot(out.2more, lower="cond-int")

out.ai <- addint(hyper, qtl=revqtl, formula=myformula)
out.ai

qtl2 <- addtoqtl(hyper, revqtl, 7, 53.6)
qtl2

qtl3 <- dropfromqtl(qtl2, index=2)
qtl3

qtl4 <- replaceqtl(hyper, qtl3, index=1, chr=1, pos=50)
qtl4

qtl5 <- reorderqtl(qtl4, c(1:3,5,4))
qtl5

stepout.a <- stepwiseqtl(hyper, additive.only=TRUE, max.qtl=6)
stepout.a

stepout.i <- stepwiseqtl(hyper, max.qtl=6)
stepout.i

############################################################
# Example 6: Internal data structure
############################################################
data(fake.bc)

class(fake.bc)

names(fake.bc)

fake.bc$pheno[1:10,]

names(fake.bc$geno)
sapply(fake.bc$geno, class)

names(fake.bc$geno[[3]])
fake.bc$geno[[3]]$data[1:5,]
fake.bc$geno[[3]]$map

names(fake.bc$geno[[3]])
fake.bc <- calc.genoprob(fake.bc, step=10, err=0.01)
names(fake.bc$geno[[3]])
fake.bc <- sim.geno(fake.bc, step=10, n.draws=8, err=0.01)
names(fake.bc$geno[[3]])
fake.bc <- argmax.geno(fake.bc, step=10, err=0.01)
names(fake.bc$geno[[3]])
fake.bc <- calc.errorlod(fake.bc, err=0.01)
names(fake.bc$geno[[3]])

names(fake.bc)
fake.bc <- est.rf(fake.bc)
names(fake.bc)

# end of rqtltour.R

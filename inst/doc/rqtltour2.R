##############################################################
# R code for "A shorter tour of R/qtl"
#
# Karl W Broman, kbroman@biostat.wisc.edu
# University of Wisconsin Madison
#
# http://www.rqtl.org
#
# 26 November 2012
##############################################################

############################################################
# Preliminaries
############################################################
library(qtl)

ls()

help(read.cross)
?read.cross

############################################################
# Data import
############################################################
sug <- read.cross("csv", "http://www.rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"), alleles=c("C", "B"))

############################################################
# Summaries
############################################################
summary(sug)

nind(sug)
nchr(sug)
totmar(sug)
nmar(sug)
nphe(sug)

plot(sug)

plotMissing(sug)
plotMap(sug)
plotPheno(sug, pheno.col=1)
plotPheno(sug, pheno.col=2)
plotPheno(sug, pheno.col=3)
plotPheno(sug, pheno.col=4)
plotPheno(sug, pheno.col=5)
plotPheno(sug, pheno.col=6)

############################################################
# Single-QTL analysis
############################################################
sug <- calc.genoprob(sug, step=1)

out.em <- scanone(sug)

summary(out.em)

summary(out.em, threshold=3)

plot(out.em)

out.hk <- scanone(sug, method="hk")

plot(out.em, out.hk, col=c("blue", "red"))

plot(out.em, col="blue")
plot(out.hk, col="red", add=TRUE)

plot(out.hk - out.em, ylim=c(-0.3, 0.3), ylab="LOD(HK)-LOD(EM)")

sug <- sim.geno(sug, step=1, n.draws=64)
out.imp <- scanone(sug, method="imp")

plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"))

plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"), chr=c(7,15))

plot(out.imp - out.em, out.hk - out.em, col=c("green", "red"), ylim=c(-1,1))

############################################################
# Permutation tests
############################################################
load(url("http://www.rqtl.org/various.RData"))

operm <- scanone(sug, method="hk", n.perm=1000)

plot(operm)

summary(operm)
summary(operm, alpha=c(0.05, 0.2))

summary(out.hk, perms=operm, alpha=0.2, pvalues=TRUE)

############################################################
# Interval estimates of QTL location
############################################################
lodint(out.hk, chr=7)
bayesint(out.hk, chr=7)

lodint(out.hk, chr=7, expandtomarkers=TRUE)
bayesint(out.hk, chr=7, expandtomarkers=TRUE)

lodint(out.hk, chr=7, drop=2)
bayesint(out.hk, chr=7, prob=0.99)

lodint(out.hk, chr=15)
bayesint(out.hk, chr=15)

############################################################
# QTL effects
############################################################
max(out.hk)
mar <- find.marker(sug, chr=7, pos=47.7)
plotPXG(sug, marker=mar)

effectplot(sug, mname1=mar)

effectplot(sug, mname1="7@47.7")

max(out.hk, chr=15)
mar2 <- find.marker(sug, chr=15, pos=12)
plotPXG(sug, marker=mar2)
effectplot(sug, mname1="15@12")

plotPXG(sug, marker=c(mar, mar2))
plotPXG(sug, marker=c(mar2, mar))

effectplot(sug, mname1="7@47.7", mname2="15@12")
effectplot(sug, mname2="7@47.7", mname1="15@12")

############################################################
# Other phenotypes
############################################################
out.hr <- scanone(sug, pheno.col=2, method="hk")

out.bw <- scanone(sug, pheno.col="bw", method="hk")

out.logbw <- scanone(sug, pheno.col=log(sug$pheno$bw), method="hk")

out.all <- scanone(sug, pheno.col=1:4, method="hk")

summary(out.all, threshold=3)

summary(out.all, threshold=3, lodcolumn=4)

summary(out.all, threshold=3, format="allpeaks")

summary(out.all, threshold=3, format="allpheno")

summary(out.all, threshold=3, format="tabByCol")
summary(out.all, threshold=3, format="tabByChr")

############################################################
# Two-dimensional, two-QTL scans
############################################################
sug <- calc.genoprob(sug, step=2)

out2 <- scantwo(sug, method="hk")

plot(out2)

plot(out2, lower="fv1")

plot(out2, lower="fv1", upper="av1")

operm2 <- scantwo(sug, method="hk", n.perm=5)

summary(out2, perms=operm2, alpha=0.2, pvalues=TRUE)

############################################################
# Multiple-QTL analyses
############################################################
sug <- calc.genoprob(sug, step=1)

qtl <- makeqtl(sug, chr=c(7,15), pos=c(47.7, 12), what="prob")

out.fq <- fitqtl(sug, qtl=qtl, method="hk")
summary(out.fq)

summary(fitqtl(sug, qtl=qtl, method="hk", get.ests=TRUE, dropone=FALSE))

out.fqi <- fitqtl(sug, qtl=qtl, method="hk", formula=y~Q1*Q2)
out.fqi <- fitqtl(sug, qtl=qtl, method="hk", formula=y~Q1+Q2+Q1:Q2)
summary(out.fqi)

addint(sug, qtl=qtl, method="hk")

rqtl <- refineqtl(sug, qtl=qtl, method="hk")
rqtl

summary(out.fqr <- fitqtl(sug, qtl=rqtl, method="hk"))

plotLodProfile(rqtl)

plot(out.hk, chr=c(7,15), col="red", add=TRUE)

out.aq <- addqtl(sug, qtl=rqtl, method="hk")

plot(out.aq)

print(pen <- calc.penalties(operm2))

out.sq <- stepwiseqtl(sug, max.qtl=5, penalties=pen, method="hk", verbose=2)
out.sq

# end of rqtltour2.R

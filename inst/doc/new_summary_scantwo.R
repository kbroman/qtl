### R code from vignette source 'new_summary_scantwo.Rnw'

###################################################
### code chunk number 1: new_summary_scantwo.Rnw:37-38
###################################################
options(width=77)


###################################################
### code chunk number 2: loaddata (eval = FALSE)
###################################################
## library(qtl)
## data(hyper)


###################################################
### code chunk number 3: loadresults
###################################################
load("hyper_results.RData")


###################################################
### code chunk number 4: scantwo (eval = FALSE)
###################################################
## hyper <- calc.genoprob(hyper, step=2.5)
## out2 <- scantwo(hyper)


###################################################
### code chunk number 5: summaryscantwoA
###################################################
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6))


###################################################
### code chunk number 6: summaryscantwoB
###################################################
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), what="full")
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), what="add")
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), what="int")


###################################################
### code chunk number 7: summaryscantwoC
###################################################
summary(out2, allpairs=FALSE)


###################################################
### code chunk number 8: summaryscantwoD
###################################################
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), df=TRUE)


###################################################
### code chunk number 9: oldsummaryscantwo
###################################################
summaryScantwoOld(out2, thresholds=c(6, 4, 4))


###################################################
### code chunk number 10: scantwoperm (eval = FALSE)
###################################################
## operm2A <- scantwo(hyper, n.perm=200)
## operm2B <- scantwo(hyper, n.perm=200)
## operm2C <- scantwo(hyper, n.perm=200)
## operm2D <- scantwo(hyper, n.perm=200)
## operm2E <- scantwo(hyper, n.perm=200)
## operm2 <- c(operm2A, operm2B, operm2C, operm2D, operm2E)


###################################################
### code chunk number 11: summaryscantwoperm
###################################################
summary(operm2, alpha=c(0.05,0.20))


###################################################
### code chunk number 12: summaryscantwopermB
###################################################
summary(out2, perms=operm2, alphas=rep(0.05, 5))


###################################################
### code chunk number 13: summaryscantwopermC
###################################################
summary(out2, perms=operm2, alphas=c(0.05, 0.05, 0, 0.05, 0.05))


###################################################
### code chunk number 14: summaryscantwoD
###################################################
summary(out2, perms=operm2, alphas=c(0.05, 0.05, 0, 0.05, 0.05), 
        pvalues=TRUE)


###################################################
### code chunk number 15: plotscantwoA (eval = FALSE)
###################################################
## plot(out2, chr=c(1,4,6,15))


###################################################
### code chunk number 16: plotscantwoAplot
###################################################
plot(out2, chr=c(1,4,6,15),layout=list(cbind(1,2),c(5,1)),
     mar1=c(4,4,0,0)+0.1, mar2=c(4,2,0,2)+0.1)


###################################################
### code chunk number 17: plotscantwoB (eval = FALSE)
###################################################
## plot(out2, chr=c(1,4,6,15), upper="cond-int")


###################################################
### code chunk number 18: plotscantwoBplot
###################################################
plot(out2, chr=c(1,4,6,15), upper="cond-int", 
     layout=list(cbind(1,2),c(5,1)),
     mar1=c(4,4,0,0)+0.1, mar2=c(4,2,0,2)+0.1)



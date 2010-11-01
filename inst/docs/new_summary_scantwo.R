###################################################
### chunk number 1: 
###################################################
#line 38 "new_summary_scantwo.Rnw"
options(width=77)


###################################################
### chunk number 2: loaddata eval=FALSE
###################################################
## #line 63 "new_summary_scantwo.Rnw"
## library(qtl)
## data(hyper)


###################################################
### chunk number 3: loadresults
###################################################
#line 68 "new_summary_scantwo.Rnw"
load("hyper_results.RData")


###################################################
### chunk number 4: scantwo eval=FALSE
###################################################
## #line 75 "new_summary_scantwo.Rnw"
## hyper <- calc.genoprob(hyper, step=2.5)
## out2 <- scantwo(hyper)


###################################################
### chunk number 5: summaryscantwoA
###################################################
#line 184 "new_summary_scantwo.Rnw"
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6))


###################################################
### chunk number 6: summaryscantwoB
###################################################
#line 207 "new_summary_scantwo.Rnw"
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), what="full")
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), what="add")
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), what="int")


###################################################
### chunk number 7: summaryscantwoC
###################################################
#line 217 "new_summary_scantwo.Rnw"
summary(out2, allpairs=FALSE)


###################################################
### chunk number 8: summaryscantwoD
###################################################
#line 224 "new_summary_scantwo.Rnw"
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), df=TRUE)


###################################################
### chunk number 9: oldsummaryscantwo
###################################################
#line 231 "new_summary_scantwo.Rnw"
summary.scantwo.old(out2, thresholds=c(6, 4, 4))


###################################################
### chunk number 10: scantwoperm eval=FALSE
###################################################
## #line 255 "new_summary_scantwo.Rnw"
## operm2A <- scantwo(hyper, n.perm=200)
## operm2B <- scantwo(hyper, n.perm=200)
## operm2C <- scantwo(hyper, n.perm=200)
## operm2D <- scantwo(hyper, n.perm=200)
## operm2E <- scantwo(hyper, n.perm=200)
## operm2 <- c(operm2A, operm2B, operm2C, operm2D, operm2E)


###################################################
### chunk number 11: summaryscantwoperm
###################################################
#line 265 "new_summary_scantwo.Rnw"
summary(operm2, alpha=c(0.05,0.20))


###################################################
### chunk number 12: summaryscantwopermB
###################################################
#line 276 "new_summary_scantwo.Rnw"
summary(out2, perms=operm2, alphas=rep(0.05, 5))


###################################################
### chunk number 13: summaryscantwopermC
###################################################
#line 282 "new_summary_scantwo.Rnw"
summary(out2, perms=operm2, alphas=c(0.05, 0.05, 0, 0.05, 0.05))


###################################################
### chunk number 14: summaryscantwoD
###################################################
#line 290 "new_summary_scantwo.Rnw"
summary(out2, perms=operm2, alphas=c(0.05, 0.05, 0, 0.05, 0.05), 
        pvalues=TRUE)


###################################################
### chunk number 15: plotscantwoA eval=FALSE
###################################################
## #line 324 "new_summary_scantwo.Rnw"
## plot(out2, chr=c(1,4,6,15))


###################################################
### chunk number 16: plotscantwoAplot
###################################################
#line 330 "new_summary_scantwo.Rnw"
plot(out2, chr=c(1,4,6,15),layout=list(cbind(1,2),c(5,1)),
     mar1=c(4,4,0,0)+0.1, mar2=c(4,2,0,2)+0.1)


###################################################
### chunk number 17: plotscantwoB eval=FALSE
###################################################
## #line 342 "new_summary_scantwo.Rnw"
## plot(out2, chr=c(1,4,6,15), upper="cond-int")


###################################################
### chunk number 18: plotscantwoBplot
###################################################
#line 348 "new_summary_scantwo.Rnw"
plot(out2, chr=c(1,4,6,15), upper="cond-int", 
     layout=list(cbind(1,2),c(5,1)),
     mar1=c(4,4,0,0)+0.1, mar2=c(4,2,0,2)+0.1)



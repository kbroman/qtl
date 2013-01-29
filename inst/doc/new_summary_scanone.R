###################################################
### chunk number 1: 
###################################################
#line 38 "new_summary_scanone.Rnw"
options(width=77)


###################################################
### chunk number 2: loaddata
###################################################
#line 66 "new_summary_scanone.Rnw"
library(qtl)
data(fake.f2)


###################################################
### chunk number 3: loadresults
###################################################
#line 71 "new_summary_scanone.Rnw"
load("fakef2_results.RData")


###################################################
### chunk number 4: scanoneA eval=FALSE
###################################################
## #line 77 "new_summary_scanone.Rnw"
## fake.f2 <- calc.genoprob(fake.f2, step=2.5)
## out.f2 <- scanone(fake.f2, method="hk")


###################################################
### chunk number 5: summaryscanoneA
###################################################
#line 85 "new_summary_scanone.Rnw"
summary(out.f2, threshold=3, df=TRUE)


###################################################
### chunk number 6: scanonepermA eval=FALSE
###################################################
## #line 103 "new_summary_scanone.Rnw"
## operm1.f2 <- scanone(fake.f2, method="hk", n.perm=500, perm.Xsp=TRUE)
## operm2.f2 <- scanone(fake.f2, method="hk", n.perm=500, perm.Xsp=TRUE)
## operm.f2 <- c(operm1.f2, operm2.f2)


###################################################
### chunk number 7: summaryscanonepermA
###################################################
#line 113 "new_summary_scanone.Rnw"
summary(operm.f2, alpha=c(0.05, 0.20))


###################################################
### chunk number 8: summaryscaononeB
###################################################
#line 120 "new_summary_scanone.Rnw"
summary(out.f2, perms=operm.f2, alpha=0.05, pvalues=TRUE)


###################################################
### chunk number 9: thestrata
###################################################
#line 131 "new_summary_scanone.Rnw"
sex <- fake.f2$pheno$sex
pgm <- fake.f2$pheno$pgm
strata <- sex + 2*pgm
table(strata)


###################################################
### chunk number 10: scanonepermB eval=FALSE
###################################################
## #line 140 "new_summary_scanone.Rnw"
## operm1.f2strat <- scanone(fake.f2, method="hk", n.perm=250,
##                          perm.Xsp=TRUE, perm.strata=strata)
## operm2.f2strat <- scanone(fake.f2, method="hk", n.perm=250,
##                          perm.Xsp=TRUE, perm.strata=strata)
## operm3.f2strat <- scanone(fake.f2, method="hk", n.perm=250,
##                          perm.Xsp=TRUE, perm.strata=strata)
## operm4.f2strat <- scanone(fake.f2, method="hk", n.perm=250,
##                          perm.Xsp=TRUE, perm.strata=strata)
## operm.f2strat <- c(operm1.f2strat, operm2.f2strat, operm3.f2strat,
##                    operm4.f2strat)


###################################################
### chunk number 11: summaryscanonepermB
###################################################
#line 154 "new_summary_scanone.Rnw"
summary(operm.f2strat, alpha=c(0.05, 0.20))


###################################################
### chunk number 12: fakebc
###################################################
#line 162 "new_summary_scanone.Rnw"
data(fake.bc)


###################################################
### chunk number 13: loaddataB
###################################################
#line 166 "new_summary_scanone.Rnw"
load("fakebc_results.RData")


###################################################
### chunk number 14: scanoneB eval=FALSE
###################################################
## #line 172 "new_summary_scanone.Rnw"
## fake.bc <- calc.genoprob(fake.bc, step=2.5)
## out.bc <- scanone(fake.bc, pheno.col=1:2, method="hk")


###################################################
### chunk number 15: summaryscanoneC
###################################################
#line 181 "new_summary_scanone.Rnw"
summary(out.bc, threshold=3)


###################################################
### chunk number 16: summaryscanoneD
###################################################
#line 187 "new_summary_scanone.Rnw"
summary(out.bc, threshold=3, lodcolumn=2)


###################################################
### chunk number 17: summaryscanoneE
###################################################
#line 193 "new_summary_scanone.Rnw"
summary(out.bc, threshold=3, format="allpheno")


###################################################
### chunk number 18: summaryscanoneF
###################################################
#line 205 "new_summary_scanone.Rnw"
summary(out.bc, threshold=c(3,2.5), format="allpeaks")


###################################################
### chunk number 19: scanonepermC eval=FALSE
###################################################
## #line 212 "new_summary_scanone.Rnw"
## operm.bc <- scanone(out.bc, pheno.col=1:2, method="hk", n.perm=1000)


###################################################
### chunk number 20: summaryscanonepermC
###################################################
#line 217 "new_summary_scanone.Rnw"
summary(operm.bc, alpha=0.05)


###################################################
### chunk number 21: summaryscanoneG
###################################################
#line 223 "new_summary_scanone.Rnw"
summary(out.bc, perms=operm.bc, alpha=0.05, format="allpeaks", 
        pvalues=TRUE)



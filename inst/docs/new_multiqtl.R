###################################################
### chunk number 1: 
###################################################
#line 40 "new_multiqtl.Rnw"
options(width=87, digits=3, scipen=4)


###################################################
### chunk number 2: myround
###################################################
#line 45 "new_multiqtl.Rnw"
source("myround.R")


###################################################
### chunk number 3: loaddata
###################################################
#line 137 "new_multiqtl.Rnw"
library(qtl)
data(hyper)


###################################################
### chunk number 4: simgeno eval=FALSE
###################################################
## #line 150 "new_multiqtl.Rnw"
## hyper <- sim.geno(hyper, step=2, n.draws=128, err=0.001)


###################################################
### chunk number 5: runsimgeno
###################################################
#line 153 "new_multiqtl.Rnw"
file <- "Rcache/simgeno.RData"
if(file.exists(file)) {
  load(file)
} else {
  set.seed(94743379)
#line 150 "new_multiqtl.Rnw#from line#158#"
hyper <- sim.geno(hyper, step=2, n.draws=128, err=0.001)
#line 159 "new_multiqtl.Rnw"
  save(hyper, file=file)
}


###################################################
### chunk number 6: makeqtl
###################################################
#line 171 "new_multiqtl.Rnw"
qtl <- makeqtl(hyper, chr=c(1, 4, 6, 15), pos=c(67.3, 30, 60, 17.5))


###################################################
### chunk number 7: print.qtl
###################################################
#line 178 "new_multiqtl.Rnw"
qtl


###################################################
### chunk number 8: plotqtl eval=FALSE
###################################################
## #line 185 "new_multiqtl.Rnw"
## plot(qtl)


###################################################
### chunk number 9: plotqtlplot
###################################################
#line 191 "new_multiqtl.Rnw"
par(mar=c(4.1,4.1,4.1,0.1), cex.axis=0.8, cex.main=1, cex=0.7)
#line 185 "new_multiqtl.Rnw#from line#192#"
plot(qtl)
#line 193 "new_multiqtl.Rnw"


###################################################
### chunk number 10: fitqtl1
###################################################
#line 208 "new_multiqtl.Rnw"
out.fq <- fitqtl(hyper, qtl=qtl, formula=y~Q1+Q2+Q3*Q4)
summary(out.fq)


###################################################
### chunk number 11: refineqtl eval=FALSE
###################################################
## #line 231 "new_multiqtl.Rnw"
## rqtl <- refineqtl(hyper, qtl=qtl, formula=y~Q1+Q2+Q3*Q4, verbose=FALSE)


###################################################
### chunk number 12: runrefineqtl
###################################################
#line 234 "new_multiqtl.Rnw"
file <- "Rcache/refineqtl.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 231 "new_multiqtl.Rnw#from line#238#"
rqtl <- refineqtl(hyper, qtl=qtl, formula=y~Q1+Q2+Q3*Q4, verbose=FALSE)
#line 239 "new_multiqtl.Rnw"
  save(rqtl, file=file)
}


###################################################
### chunk number 13: printrefineqtl
###################################################
#line 247 "new_multiqtl.Rnw"
rqtl


###################################################
### chunk number 14: fitqtl2
###################################################
#line 255 "new_multiqtl.Rnw"
out.fq2 <- fitqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4, dropone=FALSE)
summary(out.fq2)


###################################################
### chunk number 15: lodprofile eval=FALSE
###################################################
## #line 269 "new_multiqtl.Rnw"
## plotLodProfile(rqtl)


###################################################
### chunk number 16: plotlodprofile
###################################################
#line 276 "new_multiqtl.Rnw"
par(mar=c(4.1,4.1,4.1,0.1), cex=0.8)
plotLodProfile(rqtl)


###################################################
### chunk number 17: addint
###################################################
#line 318 "new_multiqtl.Rnw"
addint(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4)


###################################################
### chunk number 18: addint2
###################################################
#line 327 "new_multiqtl.Rnw"
addint(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3+Q4)


###################################################
### chunk number 19: addqtl eval=FALSE
###################################################
## #line 344 "new_multiqtl.Rnw"
## out.aq <- addqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4)


###################################################
### chunk number 20: runnaddqtl
###################################################
#line 347 "new_multiqtl.Rnw"
file <- "Rcache/addqtl.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 344 "new_multiqtl.Rnw#from line#351#"
out.aq <- addqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4)
#line 352 "new_multiqtl.Rnw"
  save(out.aq, file=file)
}


###################################################
### chunk number 21: maxaddqtl
###################################################
#line 362 "new_multiqtl.Rnw"
max(out.aq)


###################################################
### chunk number 22: plotaddqtl eval=FALSE
###################################################
## #line 368 "new_multiqtl.Rnw"
## plot(out.aq)


###################################################
### chunk number 23: plotaddqtlplot
###################################################
#line 374 "new_multiqtl.Rnw"
par(mar=c(4.1,4.1,4.1,0.1))
#line 368 "new_multiqtl.Rnw#from line#375#"
plot(out.aq)
#line 376 "new_multiqtl.Rnw"


###################################################
### chunk number 24: addqtlint eval=FALSE
###################################################
## #line 388 "new_multiqtl.Rnw"
## out.aqi <- addqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4+Q4*Q5)


###################################################
### chunk number 25: runaddqtlint
###################################################
#line 392 "new_multiqtl.Rnw"
file <- "Rcache/addqtlint.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 388 "new_multiqtl.Rnw#from line#396#"
out.aqi <- addqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4+Q4*Q5)
#line 397 "new_multiqtl.Rnw"
  save(out.aqi, file=file)
}


###################################################
### chunk number 26: plotaddqtlint eval=FALSE
###################################################
## #line 404 "new_multiqtl.Rnw"
## plot(out.aqi)


###################################################
### chunk number 27: plotaddqtlintplot
###################################################
#line 410 "new_multiqtl.Rnw"
par(mar=c(4.1,4.1,4.1,0.1))
#line 404 "new_multiqtl.Rnw#from line#411#"
plot(out.aqi)
#line 412 "new_multiqtl.Rnw"


###################################################
### chunk number 28: plotaddqtlint2 eval=FALSE
###################################################
## #line 424 "new_multiqtl.Rnw"
## plot(out.aqi - out.aq)


###################################################
### chunk number 29: plotaddqtlint2plot
###################################################
#line 430 "new_multiqtl.Rnw"
par(mar=c(4.1,4.1,4.1,0.1))
#line 424 "new_multiqtl.Rnw#from line#431#"
plot(out.aqi - out.aq)
#line 432 "new_multiqtl.Rnw"


###################################################
### chunk number 30: addpair eval=FALSE
###################################################
## #line 465 "new_multiqtl.Rnw"
## out.ap <- addpair(hyper, qtl=rqtl, chr=1, formula=y~Q2+Q3*Q4, verbose=FALSE)


###################################################
### chunk number 31: runaddpair
###################################################
#line 468 "new_multiqtl.Rnw"
file <- "Rcache/addpair.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 465 "new_multiqtl.Rnw#from line#472#"
out.ap <- addpair(hyper, qtl=rqtl, chr=1, formula=y~Q2+Q3*Q4, verbose=FALSE)
#line 473 "new_multiqtl.Rnw"
  save(out.ap, file=file)
}


###################################################
### chunk number 32: summaddpair
###################################################
#line 481 "new_multiqtl.Rnw"
summary(out.ap)


###################################################
### chunk number 33: plotaddpair eval=FALSE
###################################################
## #line 495 "new_multiqtl.Rnw"
## plot(out.ap, lower="cond-int", upper="cond-add")


###################################################
### chunk number 34: plotaddpairplot
###################################################
#line 501 "new_multiqtl.Rnw"
plot(out.ap, lower="cond-int", upper="cond-add", 
     layout=list(cbind(1,2),c(5,1)),
     mar1=c(4,4,0,0)+0.1, mar2=c(4,2,0,2)+0.1)


###################################################
### chunk number 35: addpair2 eval=FALSE
###################################################
## #line 534 "new_multiqtl.Rnw"
## out.ap2 <- addpair(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3+Q5*Q6+Q3:Q5, chr=c(7,15),
##                    verbose=FALSE)


###################################################
### chunk number 36: runaddpair2
###################################################
#line 538 "new_multiqtl.Rnw"
file <- "Rcache/addpair2.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 534 "new_multiqtl.Rnw#from line#542#"
out.ap2 <- addpair(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3+Q5*Q6+Q3:Q5, chr=c(7,15),
                   verbose=FALSE)
#line 543 "new_multiqtl.Rnw"
  save(out.ap2, file=file)
}


###################################################
### chunk number 37: summaddpair2
###################################################
#line 567 "new_multiqtl.Rnw"
summary(out.ap2)


###################################################
### chunk number 38: plotaddpair2 eval=FALSE
###################################################
## #line 592 "new_multiqtl.Rnw"
## plot(out.ap2)


###################################################
### chunk number 39: plotaddpair2plot
###################################################
#line 598 "new_multiqtl.Rnw"
plot(out.ap2, 
     layout=list(cbind(1,2),c(5,1)),
     mar1=c(4,4,0,0)+0.1, mar2=c(4,2,0,2)+0.1)


###################################################
### chunk number 40: addtoqtl
###################################################
#line 644 "new_multiqtl.Rnw"
rqtl <- addtoqtl(hyper, rqtl, 1, 43.3)
rqtl


###################################################
### chunk number 41: replaceqtl
###################################################
#line 657 "new_multiqtl.Rnw"
rqtl <- replaceqtl(hyper, rqtl, 1, 1, 73.3)
rqtl


###################################################
### chunk number 42: reorderqtl
###################################################
#line 667 "new_multiqtl.Rnw"
rqtl <- reorderqtl(rqtl, c(5,1:4))
rqtl


###################################################
### chunk number 43: dropfromqtl
###################################################
#line 674 "new_multiqtl.Rnw"
rqtl <- dropfromqtl(rqtl, 2)
rqtl


###################################################
### chunk number 44: stepqtl1 eval=FALSE
###################################################
## #line 717 "new_multiqtl.Rnw"
## stepout1 <- stepwiseqtl(hyper, additive.only=TRUE, max.qtl=6,
##                         verbose=FALSE)


###################################################
### chunk number 45: runstepqtl1
###################################################
#line 721 "new_multiqtl.Rnw"
file <- "Rcache/stepqtl1.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 717 "new_multiqtl.Rnw#from line#725#"
stepout1 <- stepwiseqtl(hyper, additive.only=TRUE, max.qtl=6,
                        verbose=FALSE)
#line 726 "new_multiqtl.Rnw"
  save(stepout1, file=file)
}


###################################################
### chunk number 46: printstepqtl1
###################################################
#line 733 "new_multiqtl.Rnw"
stepout1


###################################################
### chunk number 47: stepqtl2 eval=FALSE
###################################################
## #line 744 "new_multiqtl.Rnw"
## stepout2 <- stepwiseqtl(hyper, max.qtl=6, keeptrace=TRUE,
##                         verbose=FALSE)


###################################################
### chunk number 48: runstepqtl2
###################################################
#line 748 "new_multiqtl.Rnw"
file <- "Rcache/stepqtl2.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 744 "new_multiqtl.Rnw#from line#752#"
stepout2 <- stepwiseqtl(hyper, max.qtl=6, keeptrace=TRUE,
                        verbose=FALSE)
#line 753 "new_multiqtl.Rnw"
  save(stepout2, file=file)
}


###################################################
### chunk number 49: printstepqtl2
###################################################
#line 761 "new_multiqtl.Rnw"
stepout2


###################################################
### chunk number 50: attributenames
###################################################
#line 774 "new_multiqtl.Rnw"
names(attributes(stepout2))


###################################################
### chunk number 51: thetrace
###################################################
#line 784 "new_multiqtl.Rnw"
thetrace <- attr(stepout2, "trace")
thetrace[[1]]


###################################################
### chunk number 52: traceplot eval=FALSE
###################################################
## #line 800 "new_multiqtl.Rnw"
## par(mfrow=c(4,3))
## for(i in seq(along=thetrace))
##   plotModel(thetrace[[i]], chronly=TRUE,
##             main=paste(i, ": pLOD =", 
##               round(attr(thetrace[[i]], "pLOD"), 2)))


###################################################
### chunk number 53: plottraceplot
###################################################
#line 810 "new_multiqtl.Rnw"
par(mar=c(0.6,0.1,2.1,0.6))
#line 800 "new_multiqtl.Rnw#from line#811#"
par(mfrow=c(4,3))
for(i in seq(along=thetrace))
  plotModel(thetrace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =", 
              round(attr(thetrace[[i]], "pLOD"), 2)))
#line 812 "new_multiqtl.Rnw"



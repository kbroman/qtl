###################################################
### chunk number 1: 
###################################################
#line 39 "geneticmaps.Rnw"
options(width=87, digits=3, scipen=4)
set.seed(61777369)


###################################################
### chunk number 2: myround
###################################################
#line 46 "geneticmaps.Rnw"
source("myround.R")


###################################################
### chunk number 3: loaddata
###################################################
#line 212 "geneticmaps.Rnw"
library(qtl)
data(mapthis)


###################################################
### chunk number 4: readdata eval=FALSE
###################################################
## #line 235 "geneticmaps.Rnw"
## mapthis <- read.cross("csv", "http://www.rqtl.org/tutorials", "mapthis.csv", 
##                       estimate.map=FALSE)


###################################################
### chunk number 5: summarycross
###################################################
#line 249 "geneticmaps.Rnw"
summary(mapthis)


###################################################
### chunk number 6: plotmissing eval=FALSE
###################################################
## #line 263 "geneticmaps.Rnw"
## plot.missing(mapthis)


###################################################
### chunk number 7: plotmissingplot
###################################################
#line 269 "geneticmaps.Rnw"
par(mar=c(4.1,4.1,0.6,1.1))
plot.missing(mapthis, main="")


###################################################
### chunk number 8: plotntyped eval=FALSE
###################################################
## #line 288 "geneticmaps.Rnw"
## par(mfrow=c(1,2), las=1)
## plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
## plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", 
##      main="No. genotypes by marker")


###################################################
### chunk number 9: plotntypedplot
###################################################
#line 297 "geneticmaps.Rnw"
par(mfrow=c(1,2), las=1, cex=0.8)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", 
     main="No. genotypes by marker")


###################################################
### chunk number 10: dropind
###################################################
#line 320 "geneticmaps.Rnw"
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>50))


###################################################
### chunk number 11: dropmarkers
###################################################
#line 328 "geneticmaps.Rnw"
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 200])
mapthis <- drop.markers(mapthis, todrop)


###################################################
### chunk number 12: comparegeno eval=FALSE
###################################################
## #line 359 "geneticmaps.Rnw"
## cg <- comparegeno(mapthis)
## hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
## rug(cg[lower.tri(cg)])


###################################################
### chunk number 13: comparegenoplot
###################################################
#line 367 "geneticmaps.Rnw"
cg <- comparegeno(mapthis)
par(mar=c(4.1,4.1,0.1,0.6),las=1)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes",
     main="")
rug(cg[lower.tri(cg)])


###################################################
### chunk number 14: matchingpairs
###################################################
#line 385 "geneticmaps.Rnw"
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh


###################################################
### chunk number 15: matchinggenotypes
###################################################
#line 394 "geneticmaps.Rnw"
g <- pull.geno(mapthis)
table(g[144,], g[292,])
table(g[214,], g[216,])
table(g[238,], g[288,])


###################################################
### chunk number 16: dropmismatches
###################################################
#line 410 "geneticmaps.Rnw"
for(i in 1:nrow(wh)) {
  tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
  mapthis$geno[[1]]$data[wh[i,1],tozero] <- NA
}


###################################################
### chunk number 17: omitdup
###################################################
#line 420 "geneticmaps.Rnw"
mapthis <- subset(mapthis, ind=-wh[,2])


###################################################
### chunk number 18: finddupmar
###################################################
#line 438 "geneticmaps.Rnw"
print(dup <- findDupMarkers(mapthis, exact.only=FALSE))


###################################################
### chunk number 19: genotable
###################################################
#line 462 "geneticmaps.Rnw"
gt <- geno.table(mapthis)
gt[gt$P.value < 0.05/totmar(mapthis),]


###################################################
### chunk number 20: dropbadmarkers
###################################################
#line 475 "geneticmaps.Rnw"
todrop <- rownames(gt[gt$P.value < 1e-10,])
mapthis <- drop.markers(mapthis, todrop)


###################################################
### chunk number 21: genofreqbyind eval=FALSE
###################################################
## #line 498 "geneticmaps.Rnw"
## g <- pull.geno(mapthis)
## gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
## gfreq <- t(t(gfreq) / colSums(gfreq))
## par(mfrow=c(1,3), las=1)
## for(i in 1:3)
##   plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
##        ylim=c(0,1))


###################################################
### chunk number 22: plotgenofreqbyind
###################################################
#line 510 "geneticmaps.Rnw"
g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
       ylim=c(0,1))


###################################################
### chunk number 23: triangleplot
###################################################
#line 546 "geneticmaps.Rnw"
source("holmans_triangle.R")
par(mar=rep(0.1,4), pty="s")
triplot(labels=c("AA","AB","BB"))
tripoints(gfreq, cex=0.8)
tripoints(c(0.25, 0.5, 0.25), col="red", lwd=2, cex=1, pch=4)


###################################################
### chunk number 24: pairwiselinkage
###################################################
#line 600 "geneticmaps.Rnw"
mapthis <- est.rf(mapthis)


###################################################
### chunk number 25: checkAlleles
###################################################
#line 620 "geneticmaps.Rnw"
checkAlleles(mapthis, threshold=5)


###################################################
### chunk number 26: lodvrf eval=FALSE
###################################################
## #line 637 "geneticmaps.Rnw"
## rf <- pull.rf(mapthis)
## lod <- pull.rf(mapthis, what="lod")
## plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


###################################################
### chunk number 27: lodvrfplot
###################################################
#line 645 "geneticmaps.Rnw"
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
par(mar=c(4.1,4.1,0.6,0.6), las=1, cex=0.8)
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


###################################################
### chunk number 28: forminitialgroups
###################################################
#line 681 "geneticmaps.Rnw"
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6)
table(lg[,2])


###################################################
### chunk number 29: reorganizemarkers
###################################################
#line 703 "geneticmaps.Rnw"
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)


###################################################
### chunk number 30: plotrf eval=FALSE
###################################################
## #line 710 "geneticmaps.Rnw"
## plot.rf(mapthis, alternate.chrid=TRUE)


###################################################
### chunk number 31: plotrfplot
###################################################
#line 716 "geneticmaps.Rnw"
par(mar=c(4.1,4.1,2.1,2.1), las=1)
plot.rf(mapthis, main="", alternate.chrid=TRUE)


###################################################
### chunk number 32: plotrfonemarker eval=FALSE
###################################################
## #line 756 "geneticmaps.Rnw"
## rf <- pull.rf(mapthis)
## lod <- pull.rf(mapthis, what="lod")
## mn4 <- markernames(mapthis, chr=4)
## par(mfrow=c(2,1))
## plot(rf, mn4[3], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
## abline(h=0.5, lty=2)
## plot(lod, mn4[3], bandcol="gray70", alternate.chrid=TRUE)


###################################################
### chunk number 33: plotrfonemarkerplot
###################################################
#line 768 "geneticmaps.Rnw"
par(mar=c(4.1,4.1,1.1,0.6), las=1)
#line 756 "geneticmaps.Rnw#from line#769#"
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
mn4 <- markernames(mapthis, chr=4)
par(mfrow=c(2,1))
plot(rf, mn4[3], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, mn4[3], bandcol="gray70", alternate.chrid=TRUE)
#line 770 "geneticmaps.Rnw"


###################################################
### chunk number 34: genocrosstab
###################################################
#line 787 "geneticmaps.Rnw"
geno.crosstab(mapthis, mn4[3], mn4[1])
mn5 <- markernames(mapthis, chr=5)
geno.crosstab(mapthis, mn4[3], mn5[1])


###################################################
### chunk number 35: switchalleles
###################################################
#line 804 "geneticmaps.Rnw"
toswitch <- markernames(mapthis, chr=c(5, 7:11))
mapthis <- switchAlleles(mapthis, toswitch)


###################################################
### chunk number 36: plotrfagain eval=FALSE
###################################################
## #line 813 "geneticmaps.Rnw"
## mapthis <- est.rf(mapthis)
## plot.rf(mapthis, alternate.chrid=TRUE)


###################################################
### chunk number 37: plotrfagainplot
###################################################
#line 820 "geneticmaps.Rnw"
mapthis <- est.rf(mapthis)
par(mar=c(4.1,4.1,2.1,2.1), las=1)
plot.rf(mapthis, main="", alternate.chrid=TRUE)


###################################################
### chunk number 38: lodvrfagain eval=FALSE
###################################################
## #line 841 "geneticmaps.Rnw"
## rf <- pull.rf(mapthis)
## lod <- pull.rf(mapthis, what="lod")
## plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


###################################################
### chunk number 39: lodvrfagainplot
###################################################
#line 849 "geneticmaps.Rnw"
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
par(mar=c(4.1,4.1,0.6,0.6), las=1, cex=0.8)
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


###################################################
### chunk number 40: formgroupsagain
###################################################
#line 876 "geneticmaps.Rnw"
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6)
table(lg[,2])


###################################################
### chunk number 41: reorganizemarkersagain
###################################################
#line 883 "geneticmaps.Rnw"
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)


###################################################
### chunk number 42: plotrfyetagain eval=FALSE
###################################################
## #line 890 "geneticmaps.Rnw"
## plot.rf(mapthis)


###################################################
### chunk number 43: plotrfyetagainplot
###################################################
#line 896 "geneticmaps.Rnw"
mapthis <- est.rf(mapthis)
par(mar=c(4.1,4.1,1.6,1.6), las=1)
plot.rf(mapthis, main="")


###################################################
### chunk number 44: orderchrfive eval=FALSE
###################################################
## #line 943 "geneticmaps.Rnw"
## mapthis <- orderMarkers(mapthis, chr=5)


###################################################
### chunk number 45: orderchrfiverun
###################################################
#line 946 "geneticmaps.Rnw"
file <- "Rcache/order5.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 943 "geneticmaps.Rnw#from line#950#"
mapthis <- orderMarkers(mapthis, chr=5)
#line 951 "geneticmaps.Rnw"
  save(mapthis, file=file)
}


###################################################
### chunk number 46: chrfivemap
###################################################
#line 958 "geneticmaps.Rnw"
pull.map(mapthis, chr=5)


###################################################
### chunk number 47: ripplechr5run
###################################################
#line 992 "geneticmaps.Rnw"
file <- "Rcache/rip5.RData"
if(file.exists(file)) {
  load(file)
} else {
rip5 <- ripple(mapthis, chr=5, window=7)
save(rip5, file=file)
}


###################################################
### chunk number 48: ripplechr5 eval=FALSE
###################################################
## #line 1001 "geneticmaps.Rnw"
## rip5 <- ripple(mapthis, chr=5, window=7)


###################################################
### chunk number 49: summaryripple5
###################################################
#line 1009 "geneticmaps.Rnw"
summary(rip5)


###################################################
### chunk number 50: ripplechr5likrun
###################################################
#line 1025 "geneticmaps.Rnw"
file <- "Rcache/rip5lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip5lik <- ripple(mapthis, chr=5, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip5lik, file=file)
}


###################################################
### chunk number 51: ripplechr5lik eval=FALSE
###################################################
## #line 1035 "geneticmaps.Rnw"
## rip5lik <- ripple(mapthis, chr=5, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### chunk number 52: summaryripple5lik
###################################################
#line 1044 "geneticmaps.Rnw"
summary(rip5lik)


###################################################
### chunk number 53: compareorder
###################################################
#line 1072 "geneticmaps.Rnw"
compareorder(mapthis, chr=5, c(1:7,9,8), error.prob=0.01)
compareorder(mapthis, chr=5, c(1:7,9,8), error.prob=0.001)
compareorder(mapthis, chr=5, c(1:7,9,8), error.prob=0)


###################################################
### chunk number 54: switchorder
###################################################
#line 1088 "geneticmaps.Rnw"
mapthis <- switch.order(mapthis, chr=5, c(1:7,9,8), error.prob=0.005)
pull.map(mapthis, chr=5)


###################################################
### chunk number 55: orderchrfour eval=FALSE
###################################################
## #line 1109 "geneticmaps.Rnw"
## mapthis <- orderMarkers(mapthis, chr=4)
## pull.map(mapthis, chr=4)


###################################################
### chunk number 56: orderchrfourrun
###################################################
#line 1113 "geneticmaps.Rnw"
file <- "Rcache/order4.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 1109 "geneticmaps.Rnw#from line#1117#"
mapthis <- orderMarkers(mapthis, chr=4)
pull.map(mapthis, chr=4)
#line 1118 "geneticmaps.Rnw"
  save(mapthis, file=file)
}
pull.map(mapthis, chr=4)


###################################################
### chunk number 57: ripplechr4run
###################################################
#line 1130 "geneticmaps.Rnw"
file <- "Rcache/rip4.RData"
if(file.exists(file)) {
  load(file)
} else {
rip4 <- ripple(mapthis, chr=4, window=7)
save(rip4, file=file)
}


###################################################
### chunk number 58: ripplechr4 eval=FALSE
###################################################
## #line 1139 "geneticmaps.Rnw"
## rip4 <- ripple(mapthis, chr=4, window=7)


###################################################
### chunk number 59: summaryripple4
###################################################
#line 1147 "geneticmaps.Rnw"
summary(rip4)


###################################################
### chunk number 60: ripplechr4likrun
###################################################
#line 1155 "geneticmaps.Rnw"
file <- "Rcache/rip4lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip4lik <- ripple(mapthis, chr=4, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip4lik, file=file)
}


###################################################
### chunk number 61: ripplechr4lik eval=FALSE
###################################################
## #line 1165 "geneticmaps.Rnw"
## rip4lik <- ripple(mapthis, chr=4, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### chunk number 62: summaryripple4lik
###################################################
#line 1174 "geneticmaps.Rnw"
summary(rip4lik)


###################################################
### chunk number 63: switchmarkers4
###################################################
#line 1181 "geneticmaps.Rnw"
mapthis <- switch.order(mapthis, chr=4, c(1:8,10,9), error.prob=0.005)
pull.map(mapthis, chr=4)


###################################################
### chunk number 64: orderchrthree eval=FALSE
###################################################
## #line 1198 "geneticmaps.Rnw"
## mapthis <- orderMarkers(mapthis, chr=3)
## pull.map(mapthis, chr=3)


###################################################
### chunk number 65: orderchrthreerun
###################################################
#line 1202 "geneticmaps.Rnw"
file <- "Rcache/order3.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 1198 "geneticmaps.Rnw#from line#1206#"
mapthis <- orderMarkers(mapthis, chr=3)
pull.map(mapthis, chr=3)
#line 1207 "geneticmaps.Rnw"
  save(mapthis, file=file)
}
pull.map(mapthis, chr=3)


###################################################
### chunk number 66: ripplechr3run
###################################################
#line 1218 "geneticmaps.Rnw"
file <- "Rcache/rip3.RData"
if(file.exists(file)) {
  load(file)
} else {
rip3 <- ripple(mapthis, chr=3, window=7)
save(rip3, file=file)
}


###################################################
### chunk number 67: ripplechr3 eval=FALSE
###################################################
## #line 1227 "geneticmaps.Rnw"
## rip3 <- ripple(mapthis, chr=3, window=7)


###################################################
### chunk number 68: summaryripple3
###################################################
#line 1235 "geneticmaps.Rnw"
summary(rip3)


###################################################
### chunk number 69: ripplechr3likrun
###################################################
#line 1243 "geneticmaps.Rnw"
file <- "Rcache/rip3lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip3lik <- ripple(mapthis, chr=3, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip3lik, file=file)
}


###################################################
### chunk number 70: ripplechr3lik eval=FALSE
###################################################
## #line 1253 "geneticmaps.Rnw"
## rip3lik <- ripple(mapthis, chr=3, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### chunk number 71: summaryripple3lik
###################################################
#line 1262 "geneticmaps.Rnw"
summary(rip3lik)


###################################################
### chunk number 72: orderchrtwo eval=FALSE
###################################################
## #line 1276 "geneticmaps.Rnw"
## mapthis <- orderMarkers(mapthis, chr=2)
## pull.map(mapthis, chr=2)


###################################################
### chunk number 73: orderchrtworun
###################################################
#line 1280 "geneticmaps.Rnw"
file <- "Rcache/order2.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 1276 "geneticmaps.Rnw#from line#1284#"
mapthis <- orderMarkers(mapthis, chr=2)
pull.map(mapthis, chr=2)
#line 1285 "geneticmaps.Rnw"
  save(mapthis, file=file)
}
pull.map(mapthis, chr=2)


###################################################
### chunk number 74: ripplechr2run
###################################################
#line 1294 "geneticmaps.Rnw"
file <- "Rcache/rip2.RData"
if(file.exists(file)) {
  load(file)
} else {
rip2 <- ripple(mapthis, chr=2, window=7)
save(rip2, file=file)
}


###################################################
### chunk number 75: ripplechr2 eval=FALSE
###################################################
## #line 1303 "geneticmaps.Rnw"
## rip2 <- ripple(mapthis, chr=2, window=7)


###################################################
### chunk number 76: summaryripple2
###################################################
#line 1311 "geneticmaps.Rnw"
summary(rip2)


###################################################
### chunk number 77: ripplechr2likrun
###################################################
#line 1319 "geneticmaps.Rnw"
file <- "Rcache/rip2lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip2lik <- ripple(mapthis, chr=2, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip2lik, file=file)
}


###################################################
### chunk number 78: ripplechr2lik eval=FALSE
###################################################
## #line 1329 "geneticmaps.Rnw"
## rip2lik <- ripple(mapthis, chr=2, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### chunk number 79: summaryripple2lik
###################################################
#line 1338 "geneticmaps.Rnw"
summary(rip2lik)


###################################################
### chunk number 80: comparexo2lik eval=FALSE
###################################################
## #line 1354 "geneticmaps.Rnw"
## pat2 <- apply(rip2[,1:24], 1, paste, collapse=":")
## pat2lik <- apply(rip2lik[,1:24], 1, paste, collapse=":")
## rip2 <- rip2[match(pat2lik, pat2),]
## plot(rip2[,"obligXO"], rip2lik[,"LOD"], xlab="obligate crossover count",
##      ylab="LOD score")


###################################################
### chunk number 81: comparexo2likplot
###################################################
#line 1364 "geneticmaps.Rnw"
par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.8)
#line 1354 "geneticmaps.Rnw#from line#1365#"
pat2 <- apply(rip2[,1:24], 1, paste, collapse=":")
pat2lik <- apply(rip2lik[,1:24], 1, paste, collapse=":")
rip2 <- rip2[match(pat2lik, pat2),]
plot(rip2[,"obligXO"], rip2lik[,"LOD"], xlab="obligate crossover count",
     ylab="LOD score")
#line 1366 "geneticmaps.Rnw"


###################################################
### chunk number 82: orderchrone eval=FALSE
###################################################
## #line 1385 "geneticmaps.Rnw"
## mapthis <- orderMarkers(mapthis, chr=1)
## pull.map(mapthis, chr=1)


###################################################
### chunk number 83: orderchronerun
###################################################
#line 1389 "geneticmaps.Rnw"
file <- "Rcache/order1.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 1385 "geneticmaps.Rnw#from line#1393#"
mapthis <- orderMarkers(mapthis, chr=1)
pull.map(mapthis, chr=1)
#line 1394 "geneticmaps.Rnw"
  save(mapthis, file=file)
}
pull.map(mapthis, chr=1)


###################################################
### chunk number 84: ripplechr1run
###################################################
#line 1410 "geneticmaps.Rnw"
file <- "Rcache/rip1.RData"
if(file.exists(file)) {
  load(file)
} else {
rip1 <- ripple(mapthis, chr=1, window=7)
save(rip1, file=file)
}


###################################################
### chunk number 85: ripplechr1 eval=FALSE
###################################################
## #line 1419 "geneticmaps.Rnw"
## rip1 <- ripple(mapthis, chr=1, window=7)


###################################################
### chunk number 86: summaryripple1
###################################################
#line 1427 "geneticmaps.Rnw"
summary(rip1)


###################################################
### chunk number 87: ripplechr1likrun
###################################################
#line 1435 "geneticmaps.Rnw"
file <- "Rcache/rip1lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip1lik <- ripple(mapthis, chr=1, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip1lik, file=file)
}


###################################################
### chunk number 88: ripplechr1lik eval=FALSE
###################################################
## #line 1445 "geneticmaps.Rnw"
## rip1lik <- ripple(mapthis, chr=1, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### chunk number 89: summaryripple1lik
###################################################
#line 1454 "geneticmaps.Rnw"
summary(rip1lik)


###################################################
### chunk number 90: summarymap
###################################################
#line 1485 "geneticmaps.Rnw"
summary.map(mapthis)


###################################################
### chunk number 91: savesummarymap
###################################################
#line 1488 "geneticmaps.Rnw"
firstsummary <- summary.map(mapthis)


###################################################
### chunk number 92: plotmap eval=FALSE
###################################################
## #line 1503 "geneticmaps.Rnw"
## plot.map(mapthis, show.marker.names=TRUE)


###################################################
### chunk number 93: plotmapplot
###################################################
#line 1509 "geneticmaps.Rnw"
par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.8)
plot.map(mapthis, main="", show.marker.names=TRUE)


###################################################
### chunk number 94: plotrfonemoretime eval=FALSE
###################################################
## #line 1525 "geneticmaps.Rnw"
## plot.rf(mapthis)


###################################################
### chunk number 95: plotrfonemoretimeplot
###################################################
#line 1531 "geneticmaps.Rnw"
par(mar=c(4.1,4.1,1.6,1.6), las=1)
plot.rf(mapthis, main="")


###################################################
### chunk number 96: plotrfafterreorder eval=FALSE
###################################################
## #line 1556 "geneticmaps.Rnw"
## messedup <- switch.order(mapthis, chr=1, c(1:11,23:33,12:22), 
##                          error.prob=0.005)
## plot.rf(messedup, chr=1)


###################################################
### chunk number 97: plotrfafterreorderplot
###################################################
#line 1564 "geneticmaps.Rnw"
par(mar=c(4.1,4.1,1.6,1.6), las=1, pty="s", cex=0.8)
messedup <- switch.order(mapthis, chr=1, c(1:11,23:33,12:22), 
                         error.prob=0.005)
plot.rf(messedup, chr=1, main="")


###################################################
### chunk number 98: plotmapmessedup eval=FALSE
###################################################
## #line 1588 "geneticmaps.Rnw"
## plot.map(messedup, show.marker.names=TRUE)


###################################################
### chunk number 99: plotmapmessedupplot
###################################################
#line 1594 "geneticmaps.Rnw"
par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.8)
plot.map(messedup, main="", show.marker.names=TRUE)


###################################################
### chunk number 100: droponemarker eval=FALSE
###################################################
## #line 1636 "geneticmaps.Rnw"
## dropone <- droponemarker(mapthis, error.prob=0.005)


###################################################
### chunk number 101: droponemarkerrun
###################################################
#line 1639 "geneticmaps.Rnw"
file <- "Rcache/dropone.RData"
if(file.exists(file)) {
  load(file)
} else {
  dropone <- droponemarker(mapthis, error.prob=0.005)
  save(dropone, file=file)
}


###################################################
### chunk number 102: plotdropone eval=FALSE
###################################################
## #line 1653 "geneticmaps.Rnw"
## par(mfrow=c(2,1))
## plot(dropone, lod=1, ylim=c(-100,0))
## plot(dropone, lod=2, ylab="Change in chromosome length")


###################################################
### chunk number 103: plotdroponeplot
###################################################
#line 1661 "geneticmaps.Rnw"
par(mar=c(4.1,4.1,1.6,0.1), mfrow=c(2,1), cex=0.8)
plot(dropone, lod=1, ylim=c(-100,0))
plot(dropone, lod=2, ylab="Change in chr length (cM)")


###################################################
### chunk number 104: worstmarkers
###################################################
#line 1688 "geneticmaps.Rnw"
summary(dropone, lod.column=2)


###################################################
### chunk number 105: dropbadmarkers
###################################################
#line 1698 "geneticmaps.Rnw"
badmar <- rownames(summary(dropone, lod.column=2))[1:3]
mapthis <- drop.markers(mapthis, badmar)


###################################################
### chunk number 106: reestimatemap
###################################################
#line 1706 "geneticmaps.Rnw"
newmap <- est.map(mapthis, error.prob=0.005)
mapthis <- replace.map(mapthis, newmap)
summary.map(mapthis)


###################################################
### chunk number 107: savenewsummary
###################################################
#line 1711 "geneticmaps.Rnw"
secondsummary <- summary.map(mapthis)


###################################################
### chunk number 108: countxo eval=FALSE
###################################################
## #line 1737 "geneticmaps.Rnw"
## plot(countXO(mapthis), ylab="Number of crossovers")


###################################################
### chunk number 109: countxoplot
###################################################
#line 1743 "geneticmaps.Rnw"
par(mar=c(4.1,4.1,0.6,0.6), cex=0.8)
#line 1737 "geneticmaps.Rnw#from line#1744#"
plot(countXO(mapthis), ylab="Number of crossovers")
#line 1745 "geneticmaps.Rnw"
thecounts <- countXO(mapthis)
worst <- rev(sort(thecounts, decreasing=TRUE)[1:2])


###################################################
### chunk number 110: drophighxoind
###################################################
#line 1758 "geneticmaps.Rnw"
mapthis <- subset(mapthis, ind=(countXO(mapthis) < 50))


###################################################
### chunk number 111: rip5again
###################################################
#line 1766 "geneticmaps.Rnw"
summary(rip <- ripple(mapthis, chr=5, window=7))
summary(rip <- ripple(mapthis, chr=5, window=2, method="likelihood", 
                      error.prob=0.005))


###################################################
### chunk number 112: switchchr5again
###################################################
#line 1775 "geneticmaps.Rnw"
mapthis <- switch.order(mapthis, chr=5, c(1:7,9,8), error.prob=0.005)
pull.map(mapthis, chr=5)


###################################################
### chunk number 113: reestmapagain
###################################################
#line 1784 "geneticmaps.Rnw"
newmap <- est.map(mapthis, error.prob=0.005)
mapthis <- replace.map(mapthis, newmap)
summary.map(mapthis)


###################################################
### chunk number 114: savethirdsummary
###################################################
#line 1789 "geneticmaps.Rnw"
thirdsummary <- summary.map(mapthis)


###################################################
### chunk number 115: studyerrorrate eval=FALSE
###################################################
## #line 1816 "geneticmaps.Rnw"
## loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
## for(i in seq(along=err)) {
##   cat(i, "of", length(err), "\n")
##   tempmap <- est.map(mapthis, error.prob=err[i])
##   loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
## }
## lod <- (loglik - max(loglik))/log(10)


###################################################
### chunk number 116: runstudyerrorrate
###################################################
#line 1825 "geneticmaps.Rnw"
file <- "Rcache/errorrate.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 1816 "geneticmaps.Rnw#from line#1829#"
loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
for(i in seq(along=err)) {
  cat(i, "of", length(err), "\n")
  tempmap <- est.map(mapthis, error.prob=err[i])
  loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)
#line 1830 "geneticmaps.Rnw"
  save(err, lod, file=file)
}


###################################################
### chunk number 117: ploterrorratelik eval=FALSE
###################################################
## #line 1844 "geneticmaps.Rnw"
## plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02),
##      ylab=expression(paste(log[10], " likelihood")))


###################################################
### chunk number 118: ploterrorratelikplot
###################################################
#line 1851 "geneticmaps.Rnw"
par(mar=c(4.1,4.1,0.6,0.6), las=1)
#line 1844 "geneticmaps.Rnw#from line#1852#"
plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02),
     ylab=expression(paste(log[10], " likelihood")))
#line 1853 "geneticmaps.Rnw"


###################################################
### chunk number 119: errorlod eval=FALSE
###################################################
## #line 1889 "geneticmaps.Rnw"
## mapthis <- calc.errorlod(mapthis, error.prob=0.005)


###################################################
### chunk number 120: runerrorlod
###################################################
#line 1892 "geneticmaps.Rnw"
file <- "Rcache/errorlod.RData"
if(file.exists(file)) {
  load(file)
} else {
#line 1889 "geneticmaps.Rnw#from line#1896#"
mapthis <- calc.errorlod(mapthis, error.prob=0.005)
#line 1897 "geneticmaps.Rnw"
  save(mapthis, file=file)
}


###################################################
### chunk number 121: toperrorlod
###################################################
#line 1908 "geneticmaps.Rnw"
print(toperr <- top.errorlod(mapthis, cutoff=6))


###################################################
### chunk number 122: plotgeno eval=FALSE
###################################################
## #line 1920 "geneticmaps.Rnw"
## plot.geno(mapthis, chr=1, ind=toperr$id[toperr$chr==1],
##           cutoff=6, include.xo=FALSE)


###################################################
### chunk number 123: plotgenoplot
###################################################
#line 1927 "geneticmaps.Rnw"
par(mar=c(4.1,4.1,0.6,0.6), las=1, cex.axis=0.9)
plot.geno(mapthis, chr=1, ind=toperr$id[toperr$chr==1], main="", cex=0.8,
          include.xo=FALSE, cutoff=6)


###################################################
### chunk number 124: dropgenotypes
###################################################
#line 1953 "geneticmaps.Rnw"
mapthis.clean <- mapthis
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  mapthis.clean$geno[[chr]]$data[mapthis$pheno$id==id, mar] <- NA
}


###################################################
### chunk number 125: segdis eval=FALSE
###################################################
## #line 1972 "geneticmaps.Rnw"
## gt <- geno.table(mapthis, scanone.output=TRUE)
## par(mfrow=c(2,1))
## plot(gt, ylab=expression(paste(-log[10], " P-value")))
## plot(gt, lod=3:5, ylab="Genotype frequency")
## abline(h=c(0.25, 0.5), lty=2, col="gray")


###################################################
### chunk number 126: plotsegdis
###################################################
#line 1982 "geneticmaps.Rnw"
gt <- geno.table(mapthis, scanone.output=TRUE)
par(mar=c(4.1,4.1,0.6,0.6), las=1, mfrow=c(2,1), cex=0.8)
plot(gt, ylab=expression(paste(-log[10], " P-value")))
plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")


###################################################
### chunk number 127: plotfinalmap eval=FALSE
###################################################
## #line 2026 "geneticmaps.Rnw"
## plot.map(mapthis, show.marker.names=TRUE)


###################################################
### chunk number 128: plotfinalmapplot
###################################################
#line 2032 "geneticmaps.Rnw"
par(las=1, mar=c(4.6,4.6,0.6,0.6), cex=0.8)
plot.map(mapthis, main="", show.marker.names=TRUE)



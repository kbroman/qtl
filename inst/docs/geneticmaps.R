### R code from vignette source 'geneticmaps.Rnw'

###################################################
### code chunk number 1: geneticmaps.Rnw:38-40
###################################################
options(width=87, digits=3, scipen=4)
set.seed(61777369)


###################################################
### code chunk number 2: myround
###################################################
source("myround.R")


###################################################
### code chunk number 3: loaddata
###################################################
library(qtl)
data(mapthis)


###################################################
### code chunk number 4: readdata (eval = FALSE)
###################################################
## mapthis <- read.cross("csv", "http://www.rqtl.org/tutorials", "mapthis.csv", 
##                       estimate.map=FALSE)


###################################################
### code chunk number 5: summarycross
###################################################
summary(mapthis)


###################################################
### code chunk number 6: plotmissing (eval = FALSE)
###################################################
## plotMissing(mapthis)


###################################################
### code chunk number 7: plotmissingplot
###################################################
par(mar=c(4.1,4.1,0.6,1.1))
plotMissing(mapthis, main="")


###################################################
### code chunk number 8: plotntyped (eval = FALSE)
###################################################
## par(mfrow=c(1,2), las=1)
## plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
## plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", 
##      main="No. genotypes by marker")


###################################################
### code chunk number 9: plotntypedplot
###################################################
par(mfrow=c(1,2), las=1, cex=0.8)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", 
     main="No. genotypes by marker")


###################################################
### code chunk number 10: dropind
###################################################
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>50))


###################################################
### code chunk number 11: dropmarkers
###################################################
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 200])
mapthis <- drop.markers(mapthis, todrop)


###################################################
### code chunk number 12: comparegeno (eval = FALSE)
###################################################
## cg <- comparegeno(mapthis)
## hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
## rug(cg[lower.tri(cg)])


###################################################
### code chunk number 13: comparegenoplot
###################################################
cg <- comparegeno(mapthis)
par(mar=c(4.1,4.1,0.1,0.6),las=1)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes",
     main="")
rug(cg[lower.tri(cg)])


###################################################
### code chunk number 14: matchingpairs
###################################################
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh


###################################################
### code chunk number 15: matchinggenotypes
###################################################
g <- pull.geno(mapthis)
table(g[144,], g[292,])
table(g[214,], g[216,])
table(g[238,], g[288,])


###################################################
### code chunk number 16: dropmismatches
###################################################
for(i in 1:nrow(wh)) {
  tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
  mapthis$geno[[1]]$data[wh[i,1],tozero] <- NA
}


###################################################
### code chunk number 17: omitdup
###################################################
mapthis <- subset(mapthis, ind=-wh[,2])


###################################################
### code chunk number 18: finddupmar
###################################################
print(dup <- findDupMarkers(mapthis, exact.only=FALSE))


###################################################
### code chunk number 19: genotable
###################################################
gt <- geno.table(mapthis)
gt[gt$P.value < 0.05/totmar(mapthis),]


###################################################
### code chunk number 20: dropbadmarkers
###################################################
todrop <- rownames(gt[gt$P.value < 1e-10,])
mapthis <- drop.markers(mapthis, todrop)


###################################################
### code chunk number 21: genofreqbyind (eval = FALSE)
###################################################
## g <- pull.geno(mapthis)
## gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
## gfreq <- t(t(gfreq) / colSums(gfreq))
## par(mfrow=c(1,3), las=1)
## for(i in 1:3)
##   plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
##        ylim=c(0,1))


###################################################
### code chunk number 22: plotgenofreqbyind
###################################################
g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
       ylim=c(0,1))


###################################################
### code chunk number 23: triangleplot
###################################################
source("holmans_triangle.R")
par(mar=rep(0.1,4), pty="s")
triplot(labels=c("AA","AB","BB"))
tripoints(gfreq, cex=0.8)
tripoints(c(0.25, 0.5, 0.25), col="red", lwd=2, cex=1, pch=4)


###################################################
### code chunk number 24: pairwiselinkage
###################################################
mapthis <- est.rf(mapthis)


###################################################
### code chunk number 25: checkAlleles
###################################################
checkAlleles(mapthis, threshold=5)


###################################################
### code chunk number 26: lodvrf (eval = FALSE)
###################################################
## rf <- pull.rf(mapthis)
## lod <- pull.rf(mapthis, what="lod")
## plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


###################################################
### code chunk number 27: lodvrfplot
###################################################
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
par(mar=c(4.1,4.1,0.6,0.6), las=1, cex=0.8)
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


###################################################
### code chunk number 28: forminitialgroups
###################################################
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6)
table(lg[,2])


###################################################
### code chunk number 29: reorganizemarkers
###################################################
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)


###################################################
### code chunk number 30: plotrf (eval = FALSE)
###################################################
## plotRF(mapthis, alternate.chrid=TRUE)


###################################################
### code chunk number 31: plotrfplot
###################################################
par(mar=c(4.1,4.1,2.1,2.1), las=1)
plotRF(mapthis, main="", alternate.chrid=TRUE)


###################################################
### code chunk number 32: plotrfonemarker (eval = FALSE)
###################################################
## rf <- pull.rf(mapthis)
## lod <- pull.rf(mapthis, what="lod")
## mn4 <- markernames(mapthis, chr=4)
## par(mfrow=c(2,1))
## plot(rf, mn4[3], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
## abline(h=0.5, lty=2)
## plot(lod, mn4[3], bandcol="gray70", alternate.chrid=TRUE)


###################################################
### code chunk number 33: plotrfonemarkerplot
###################################################
par(mar=c(4.1,4.1,1.1,0.6), las=1)
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
mn4 <- markernames(mapthis, chr=4)
par(mfrow=c(2,1))
plot(rf, mn4[3], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, mn4[3], bandcol="gray70", alternate.chrid=TRUE)


###################################################
### code chunk number 34: genocrosstab
###################################################
geno.crosstab(mapthis, mn4[3], mn4[1])
mn5 <- markernames(mapthis, chr=5)
geno.crosstab(mapthis, mn4[3], mn5[1])


###################################################
### code chunk number 35: switchalleles
###################################################
toswitch <- markernames(mapthis, chr=c(5, 7:11))
mapthis <- switchAlleles(mapthis, toswitch)


###################################################
### code chunk number 36: plotrfagain (eval = FALSE)
###################################################
## mapthis <- est.rf(mapthis)
## plotRF(mapthis, alternate.chrid=TRUE)


###################################################
### code chunk number 37: plotrfagainplot
###################################################
mapthis <- est.rf(mapthis)
par(mar=c(4.1,4.1,2.1,2.1), las=1)
plotRF(mapthis, main="", alternate.chrid=TRUE)


###################################################
### code chunk number 38: lodvrfagain (eval = FALSE)
###################################################
## rf <- pull.rf(mapthis)
## lod <- pull.rf(mapthis, what="lod")
## plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


###################################################
### code chunk number 39: lodvrfagainplot
###################################################
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
par(mar=c(4.1,4.1,0.6,0.6), las=1, cex=0.8)
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


###################################################
### code chunk number 40: formgroupsagain
###################################################
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6)
table(lg[,2])


###################################################
### code chunk number 41: reorganizemarkersagain
###################################################
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)


###################################################
### code chunk number 42: plotrfyetagain (eval = FALSE)
###################################################
## plotRF(mapthis)


###################################################
### code chunk number 43: plotrfyetagainplot
###################################################
mapthis <- est.rf(mapthis)
par(mar=c(4.1,4.1,1.6,1.6), las=1)
plotRF(mapthis, main="")


###################################################
### code chunk number 44: orderchrfive (eval = FALSE)
###################################################
## mapthis <- orderMarkers(mapthis, chr=5)


###################################################
### code chunk number 45: orderchrfiverun
###################################################
file <- "Rcache/order5.RData"
if(file.exists(file)) {
  load(file)
} else {
mapthis <- orderMarkers(mapthis, chr=5)
  save(mapthis, file=file)
}


###################################################
### code chunk number 46: chrfivemap
###################################################
pull.map(mapthis, chr=5)


###################################################
### code chunk number 47: ripplechr5run
###################################################
file <- "Rcache/rip5.RData"
if(file.exists(file)) {
  load(file)
} else {
rip5 <- ripple(mapthis, chr=5, window=7)
save(rip5, file=file)
}


###################################################
### code chunk number 48: ripplechr5 (eval = FALSE)
###################################################
## rip5 <- ripple(mapthis, chr=5, window=7)


###################################################
### code chunk number 49: summaryripple5
###################################################
summary(rip5)


###################################################
### code chunk number 50: ripplechr5likrun
###################################################
file <- "Rcache/rip5lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip5lik <- ripple(mapthis, chr=5, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip5lik, file=file)
}


###################################################
### code chunk number 51: ripplechr5lik (eval = FALSE)
###################################################
## rip5lik <- ripple(mapthis, chr=5, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### code chunk number 52: summaryripple5lik
###################################################
summary(rip5lik)


###################################################
### code chunk number 53: compareorder
###################################################
compareorder(mapthis, chr=5, c(1:7,9,8), error.prob=0.01)
compareorder(mapthis, chr=5, c(1:7,9,8), error.prob=0.001)
compareorder(mapthis, chr=5, c(1:7,9,8), error.prob=0)


###################################################
### code chunk number 54: switchorder
###################################################
mapthis <- switch.order(mapthis, chr=5, c(1:7,9,8), error.prob=0.005)
pull.map(mapthis, chr=5)


###################################################
### code chunk number 55: orderchrfour (eval = FALSE)
###################################################
## mapthis <- orderMarkers(mapthis, chr=4)
## pull.map(mapthis, chr=4)


###################################################
### code chunk number 56: orderchrfourrun
###################################################
file <- "Rcache/order4.RData"
if(file.exists(file)) {
  load(file)
} else {
mapthis <- orderMarkers(mapthis, chr=4)
pull.map(mapthis, chr=4)
  save(mapthis, file=file)
}
pull.map(mapthis, chr=4)


###################################################
### code chunk number 57: ripplechr4run
###################################################
file <- "Rcache/rip4.RData"
if(file.exists(file)) {
  load(file)
} else {
rip4 <- ripple(mapthis, chr=4, window=7)
save(rip4, file=file)
}


###################################################
### code chunk number 58: ripplechr4 (eval = FALSE)
###################################################
## rip4 <- ripple(mapthis, chr=4, window=7)


###################################################
### code chunk number 59: summaryripple4
###################################################
summary(rip4)


###################################################
### code chunk number 60: ripplechr4likrun
###################################################
file <- "Rcache/rip4lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip4lik <- ripple(mapthis, chr=4, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip4lik, file=file)
}


###################################################
### code chunk number 61: ripplechr4lik (eval = FALSE)
###################################################
## rip4lik <- ripple(mapthis, chr=4, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### code chunk number 62: summaryripple4lik
###################################################
summary(rip4lik)


###################################################
### code chunk number 63: switchmarkers4
###################################################
mapthis <- switch.order(mapthis, chr=4, c(1:8,10,9), error.prob=0.005)
pull.map(mapthis, chr=4)


###################################################
### code chunk number 64: orderchrthree (eval = FALSE)
###################################################
## mapthis <- orderMarkers(mapthis, chr=3)
## pull.map(mapthis, chr=3)


###################################################
### code chunk number 65: orderchrthreerun
###################################################
file <- "Rcache/order3.RData"
if(file.exists(file)) {
  load(file)
} else {
mapthis <- orderMarkers(mapthis, chr=3)
pull.map(mapthis, chr=3)
  save(mapthis, file=file)
}
pull.map(mapthis, chr=3)


###################################################
### code chunk number 66: ripplechr3run
###################################################
file <- "Rcache/rip3.RData"
if(file.exists(file)) {
  load(file)
} else {
rip3 <- ripple(mapthis, chr=3, window=7)
save(rip3, file=file)
}


###################################################
### code chunk number 67: ripplechr3 (eval = FALSE)
###################################################
## rip3 <- ripple(mapthis, chr=3, window=7)


###################################################
### code chunk number 68: summaryripple3
###################################################
summary(rip3)


###################################################
### code chunk number 69: ripplechr3likrun
###################################################
file <- "Rcache/rip3lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip3lik <- ripple(mapthis, chr=3, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip3lik, file=file)
}


###################################################
### code chunk number 70: ripplechr3lik (eval = FALSE)
###################################################
## rip3lik <- ripple(mapthis, chr=3, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### code chunk number 71: summaryripple3lik
###################################################
summary(rip3lik)


###################################################
### code chunk number 72: orderchrtwo (eval = FALSE)
###################################################
## mapthis <- orderMarkers(mapthis, chr=2)
## pull.map(mapthis, chr=2)


###################################################
### code chunk number 73: orderchrtworun
###################################################
file <- "Rcache/order2.RData"
if(file.exists(file)) {
  load(file)
} else {
mapthis <- orderMarkers(mapthis, chr=2)
pull.map(mapthis, chr=2)
  save(mapthis, file=file)
}
pull.map(mapthis, chr=2)


###################################################
### code chunk number 74: ripplechr2run
###################################################
file <- "Rcache/rip2.RData"
if(file.exists(file)) {
  load(file)
} else {
rip2 <- ripple(mapthis, chr=2, window=7)
save(rip2, file=file)
}


###################################################
### code chunk number 75: ripplechr2 (eval = FALSE)
###################################################
## rip2 <- ripple(mapthis, chr=2, window=7)


###################################################
### code chunk number 76: summaryripple2
###################################################
summary(rip2)


###################################################
### code chunk number 77: ripplechr2likrun
###################################################
file <- "Rcache/rip2lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip2lik <- ripple(mapthis, chr=2, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip2lik, file=file)
}


###################################################
### code chunk number 78: ripplechr2lik (eval = FALSE)
###################################################
## rip2lik <- ripple(mapthis, chr=2, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### code chunk number 79: summaryripple2lik
###################################################
summary(rip2lik)


###################################################
### code chunk number 80: comparexo2lik (eval = FALSE)
###################################################
## pat2 <- apply(rip2[,1:24], 1, paste, collapse=":")
## pat2lik <- apply(rip2lik[,1:24], 1, paste, collapse=":")
## rip2 <- rip2[match(pat2lik, pat2),]
## plot(rip2[,"obligXO"], rip2lik[,"LOD"], xlab="obligate crossover count",
##      ylab="LOD score")


###################################################
### code chunk number 81: comparexo2likplot
###################################################
par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.8)
pat2 <- apply(rip2[,1:24], 1, paste, collapse=":")
pat2lik <- apply(rip2lik[,1:24], 1, paste, collapse=":")
rip2 <- rip2[match(pat2lik, pat2),]
plot(rip2[,"obligXO"], rip2lik[,"LOD"], xlab="obligate crossover count",
     ylab="LOD score")


###################################################
### code chunk number 82: orderchrone (eval = FALSE)
###################################################
## mapthis <- orderMarkers(mapthis, chr=1)
## pull.map(mapthis, chr=1)


###################################################
### code chunk number 83: orderchronerun
###################################################
file <- "Rcache/order1.RData"
if(file.exists(file)) {
  load(file)
} else {
mapthis <- orderMarkers(mapthis, chr=1)
pull.map(mapthis, chr=1)
  save(mapthis, file=file)
}
pull.map(mapthis, chr=1)


###################################################
### code chunk number 84: ripplechr1run
###################################################
file <- "Rcache/rip1.RData"
if(file.exists(file)) {
  load(file)
} else {
rip1 <- ripple(mapthis, chr=1, window=7)
save(rip1, file=file)
}


###################################################
### code chunk number 85: ripplechr1 (eval = FALSE)
###################################################
## rip1 <- ripple(mapthis, chr=1, window=7)


###################################################
### code chunk number 86: summaryripple1
###################################################
summary(rip1)


###################################################
### code chunk number 87: ripplechr1likrun
###################################################
file <- "Rcache/rip1lik.RData"
if(file.exists(file)) {
  load(file)
} else {
rip1lik <- ripple(mapthis, chr=1, window=4, method="likelihood", 
                  error.prob=0.005)
save(rip1lik, file=file)
}


###################################################
### code chunk number 88: ripplechr1lik (eval = FALSE)
###################################################
## rip1lik <- ripple(mapthis, chr=1, window=4, method="likelihood", 
##                   error.prob=0.005)


###################################################
### code chunk number 89: summaryripple1lik
###################################################
summary(rip1lik)


###################################################
### code chunk number 90: summarymap
###################################################
summaryMap(mapthis)


###################################################
### code chunk number 91: savesummarymap
###################################################
firstsummary <- summaryMap(mapthis)


###################################################
### code chunk number 92: plotmap (eval = FALSE)
###################################################
## plotMap(mapthis, show.marker.names=TRUE)


###################################################
### code chunk number 93: plotmapplot
###################################################
par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.8)
plotMap(mapthis, main="", show.marker.names=TRUE)


###################################################
### code chunk number 94: plotrfonemoretime (eval = FALSE)
###################################################
## plotRF(mapthis)


###################################################
### code chunk number 95: plotrfonemoretimeplot
###################################################
par(mar=c(4.1,4.1,1.6,1.6), las=1)
plotRF(mapthis, main="")


###################################################
### code chunk number 96: plotrfafterreorder (eval = FALSE)
###################################################
## messedup <- switch.order(mapthis, chr=1, c(1:11,23:33,12:22), 
##                          error.prob=0.005)
## plotRF(messedup, chr=1)


###################################################
### code chunk number 97: plotrfafterreorderplot
###################################################
par(mar=c(4.1,4.1,1.6,1.6), las=1, pty="s", cex=0.8)
messedup <- switch.order(mapthis, chr=1, c(1:11,23:33,12:22), 
                         error.prob=0.005)
plotRF(messedup, chr=1, main="")


###################################################
### code chunk number 98: plotmapmessedup (eval = FALSE)
###################################################
## plotMap(messedup, show.marker.names=TRUE)


###################################################
### code chunk number 99: plotmapmessedupplot
###################################################
par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.8)
plotMap(messedup, main="", show.marker.names=TRUE)


###################################################
### code chunk number 100: droponemarker (eval = FALSE)
###################################################
## dropone <- droponemarker(mapthis, error.prob=0.005)


###################################################
### code chunk number 101: droponemarkerrun
###################################################
file <- "Rcache/dropone.RData"
if(file.exists(file)) {
  load(file)
} else {
  dropone <- droponemarker(mapthis, error.prob=0.005)
  save(dropone, file=file)
}


###################################################
### code chunk number 102: plotdropone (eval = FALSE)
###################################################
## par(mfrow=c(2,1))
## plot(dropone, lod=1, ylim=c(-100,0))
## plot(dropone, lod=2, ylab="Change in chromosome length")


###################################################
### code chunk number 103: plotdroponeplot
###################################################
par(mar=c(4.1,4.1,1.6,0.1), mfrow=c(2,1), cex=0.8)
plot(dropone, lod=1, ylim=c(-100,0))
plot(dropone, lod=2, ylab="Change in chr length (cM)")


###################################################
### code chunk number 104: worstmarkers
###################################################
summary(dropone, lod.column=2)


###################################################
### code chunk number 105: dropbadmarkers
###################################################
badmar <- rownames(summary(dropone, lod.column=2))[1:3]
mapthis <- drop.markers(mapthis, badmar)


###################################################
### code chunk number 106: reestimatemap
###################################################
newmap <- est.map(mapthis, error.prob=0.005)
mapthis <- replace.map(mapthis, newmap)
summaryMap(mapthis)


###################################################
### code chunk number 107: savenewsummary
###################################################
secondsummary <- summaryMap(mapthis)


###################################################
### code chunk number 108: countxo (eval = FALSE)
###################################################
## plot(countXO(mapthis), ylab="Number of crossovers")


###################################################
### code chunk number 109: countxoplot
###################################################
par(mar=c(4.1,4.1,0.6,0.6), cex=0.8)
plot(countXO(mapthis), ylab="Number of crossovers")
thecounts <- countXO(mapthis)
worst <- rev(sort(thecounts, decreasing=TRUE)[1:2])


###################################################
### code chunk number 110: drophighxoind
###################################################
mapthis <- subset(mapthis, ind=(countXO(mapthis) < 50))


###################################################
### code chunk number 111: rip5again
###################################################
summary(rip <- ripple(mapthis, chr=5, window=7))
summary(rip <- ripple(mapthis, chr=5, window=2, method="likelihood", 
                      error.prob=0.005))


###################################################
### code chunk number 112: switchchr5again
###################################################
mapthis <- switch.order(mapthis, chr=5, c(1:7,9,8), error.prob=0.005)
pull.map(mapthis, chr=5)


###################################################
### code chunk number 113: reestmapagain
###################################################
newmap <- est.map(mapthis, error.prob=0.005)
mapthis <- replace.map(mapthis, newmap)
summaryMap(mapthis)


###################################################
### code chunk number 114: savethirdsummary
###################################################
thirdsummary <- summaryMap(mapthis)


###################################################
### code chunk number 115: studyerrorrate (eval = FALSE)
###################################################
## loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
## for(i in seq(along=err)) {
##   cat(i, "of", length(err), "\n")
##   tempmap <- est.map(mapthis, error.prob=err[i])
##   loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
## }
## lod <- (loglik - max(loglik))/log(10)


###################################################
### code chunk number 116: runstudyerrorrate
###################################################
file <- "Rcache/errorrate.RData"
if(file.exists(file)) {
  load(file)
} else {
loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
for(i in seq(along=err)) {
  cat(i, "of", length(err), "\n")
  tempmap <- est.map(mapthis, error.prob=err[i])
  loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)
  save(err, lod, file=file)
}


###################################################
### code chunk number 117: ploterrorratelik (eval = FALSE)
###################################################
## plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02),
##      ylab=expression(paste(log[10], " likelihood")))


###################################################
### code chunk number 118: ploterrorratelikplot
###################################################
par(mar=c(4.1,4.1,0.6,0.6), las=1)
plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02),
     ylab=expression(paste(log[10], " likelihood")))


###################################################
### code chunk number 119: errorlod (eval = FALSE)
###################################################
## mapthis <- calc.errorlod(mapthis, error.prob=0.005)


###################################################
### code chunk number 120: runerrorlod
###################################################
file <- "Rcache/errorlod.RData"
if(file.exists(file)) {
  load(file)
} else {
mapthis <- calc.errorlod(mapthis, error.prob=0.005)
  save(mapthis, file=file)
}


###################################################
### code chunk number 121: toperrorlod
###################################################
print(toperr <- top.errorlod(mapthis, cutoff=6))


###################################################
### code chunk number 122: plotgeno (eval = FALSE)
###################################################
## plotGeno(mapthis, chr=1, ind=toperr$id[toperr$chr==1],
##           cutoff=6, include.xo=FALSE)


###################################################
### code chunk number 123: plotgenoplot
###################################################
par(mar=c(4.1,4.1,0.6,0.6), las=1, cex.axis=0.9)
plotGeno(mapthis, chr=1, ind=toperr$id[toperr$chr==1], main="", cex=0.8,
          include.xo=FALSE, cutoff=6)


###################################################
### code chunk number 124: dropgenotypes
###################################################
mapthis.clean <- mapthis
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  mapthis.clean$geno[[chr]]$data[mapthis$pheno$id==id, mar] <- NA
}


###################################################
### code chunk number 125: segdis (eval = FALSE)
###################################################
## gt <- geno.table(mapthis, scanone.output=TRUE)
## par(mfrow=c(2,1))
## plot(gt, ylab=expression(paste(-log[10], " P-value")))
## plot(gt, lod=3:5, ylab="Genotype frequency")
## abline(h=c(0.25, 0.5), lty=2, col="gray")


###################################################
### code chunk number 126: plotsegdis
###################################################
gt <- geno.table(mapthis, scanone.output=TRUE)
par(mar=c(4.1,4.1,0.6,0.6), las=1, mfrow=c(2,1), cex=0.8)
plot(gt, ylab=expression(paste(-log[10], " P-value")))
plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")


###################################################
### code chunk number 127: plotfinalmap (eval = FALSE)
###################################################
## plotMap(mapthis, show.marker.names=TRUE)


###################################################
### code chunk number 128: plotfinalmapplot
###################################################
par(las=1, mar=c(4.6,4.6,0.6,0.6), cex=0.8)
plotMap(mapthis, main="", show.marker.names=TRUE)



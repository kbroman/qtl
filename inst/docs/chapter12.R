###################################################
### chunk number 1: 
###################################################
library(qtl)
data(map10)
mycross <- sim.cross(map10, type="f2", n.ind=100, missing.prob=0.02)
plot.geno(mycross)


###################################################
### chunk number 2: 
###################################################
#warns because mqm doesnt handle the X chromosome yet
augmentedcross <- mqmaugment(mycross, augment_aboveprob=1)
plot.geno(augmentedcross)


###################################################
### chunk number 3: 
###################################################
augmentedcross <- mqmaugment(mycross, augment_aboveprob=10)
plot.geno(augmentedcross)


###################################################
### chunk number 4: 
###################################################
data(hyper)
colors <- c("Black","Green")
lines <- c(2,1)
h_no_missing <- mqmaugment(hyper, augment_aboveprob=1)
h_no_missing <- calc.errorlod()
result <- mqm(h_no_missing)
result_compare <- scanone(h_no_missing)


###################################################
### chunk number 5: 
###################################################
plot(result, result_compare, col=colors, lwd=lines)


###################################################
### chunk number 6: 
###################################################
#Summary results shows the highest lodscores per chromosome
summary(result)
#find.marker extracts the marker nearest to chr4 at 30Cm
find.marker(h_no_missing,chr=4,pos=30)
#which.marker translates the name into a cofacotr number
toset <- which.marker(h_no_missing,"D4Mit164")
cofactorlist <- mqmcofactors(h_no_missing,toset)
#scan again
result <- mqm(h_no_missing, cofactorlist,plot=T)


###################################################
### chunk number 7: 
###################################################
plot(result, result_compare, col=colors, lwd=lines)


###################################################
### chunk number 8: 
###################################################
summary(result)
#we can combine find marker and which.marker commands and expand our toset variable
toset <- c(toset,which.marker(h_no_missing,find.marker(h_no_missing,1,70)))
cofactorlist <- mqmcofactors(h_no_missing,toset)
#scan again
result <- mqm(h_no_missing, cofactorlist,plot=T)


###################################################
### chunk number 9: 
###################################################
plot(result, result_compare, col=colors, lwd=lines)


###################################################
### chunk number 10: 
###################################################
cofactorlist <- mqmcofactorsEach(h_no_missing,5)
result <- mqm(h_no_missing, cofactorlist, plot=T)


###################################################
### chunk number 11: 
###################################################
plot(result, result_compare, col=colors, lwd=lines)


###################################################
### chunk number 12: 
###################################################
result <- mqm(h_no_missing, cofactorlist, alfa=0.002, plot=T)


###################################################
### chunk number 13: 
###################################################
plot.pxg(h_no_missing,marker="D1Mit102")


###################################################
### chunk number 14: 
###################################################
effectplot(h_no_missing, mname1="D1Mit19", mname2="D1Mit102")


###################################################
### chunk number 15: 
###################################################
effectplot(h_no_missing, mname1="D1Mit102", mname2="D5Mit213")


###################################################
### chunk number 16: 
###################################################
library(snow)
results <- bootstrap(h_no_missing,mqm,cofactors=cofactorlist,plot=T,verbose=T,n.clusters=2,n.run=25,b.size=25)


###################################################
### chunk number 17: 
###################################################
data(multitrait)
multifilled <- fill.geno(multitrait)
resall <- mqmall(multifilled,n.clusters=2,verbose=T)
mqmplotall(resall,"I")


###################################################
### chunk number 18: 
###################################################
cofactorlist <- mqmcofactorsEach(multifilled,3)
resall <- mqmall(multifilled,cofactors=cofactorlist,n.clusters=2,verbose=T)
mqmplotall(resall,"I")


###################################################
### chunk number 19: 
###################################################
mqmplotnice(resall,legendloc=1)



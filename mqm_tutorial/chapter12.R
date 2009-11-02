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
result <- mqm(h_no_missing)
result_compare <- scanone(h_no_missing)


###################################################
### chunk number 5: 
###################################################
plot(result, result_compare, col=colors, lwd=lines)


###################################################
### chunk number 6: 
###################################################
summary(result)
find.marker(h_no_missing,4,30)
toset <- which.marker(h_no_missing,"D4Mit164")
cofactorlist <- mqmcofactors(h_no_missing,toset)
#scan again
result <- mqm(h_no_missing, cofactorlist)
plot(result, result_compare, col=colors, lwd=lines)


###################################################
### chunk number 7: 
###################################################
summary(result)
#we can combine find marker and which.marker commands and expand our toset variable
toset <- c(toset,which.marker(h_no_missing,find.marker(h_no_missing,1,70)))
cofactorlist <- mqmcofactors(h_no_missing,toset)
#scan again
result <- mqm(h_no_missing, cofactorlist)


###################################################
### chunk number 8: 
###################################################
plot(result, result_compare, col=colors, lwd=lines)


###################################################
### chunk number 9: 
###################################################
cofactorlist <- mqmcofactorsEach(h_no_missing,5)
result <- mqm(h_no_missing, cofactorlist, plot=T)


###################################################
### chunk number 10: 
###################################################
#We could plot also back against the scanone function.
plot(result, result_compare, col=colors, lwd=lines)


###################################################
### chunk number 11: 
###################################################
result <- mqm(h_no_missing, cofactorlist, alfa=0.002, plot=T)


###################################################
### chunk number 12: 
###################################################
plot.pxg(h_no_missing,marker="D1Mit102")


###################################################
### chunk number 13: 
###################################################
effectplot(h_no_missing, mname1="D1Mit19", mname2="D1Mit102")


###################################################
### chunk number 14: 
###################################################
effectplot(h_no_missing, mname1="D1Mit102", mname2="D5Mit213")


###################################################
### chunk number 15: 
###################################################
library(snow)
results <- bootstrap(h_no_missing,mqm,cofactors=cofactorlist,plot=T,n.run=2,b.size=1)



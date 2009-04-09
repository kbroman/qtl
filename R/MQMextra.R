extra <- function(a,b,x){
	.C("R_gammaln",as.double(3.0))
	lgamma(a)

	.C("R_betacf",as.double(3.0),as.double(3.0),as.double(0.9))

	.C("r_betai",a= as.double(3.0),b= as.double(3.0),x= as.double(0.9))
	pbeta(x,a,b)
	
	.C("R_Lnorm",as.double(3),as.double(8))
	
	resm <- NULL
	resm2 <- NULL
	resm3 <- NULL
	for(x in seq(-1,1,0.01)){
		res <- NULL
		res2 <- NULL
		res3 <- NULL
		for(y in 1:10){
			aa <- .C("R_Lnorm",as.double(x),as.double(y))
			bb <- dnorm(x,0,y)
			cc <- Lnormal(x,y)
			res <- c(res,aa[2][[1]])
			res2 <- c(res2,bb)
			res3 <- c(res3,cc)
		}
		resm <- rbind(resm,res)
		resm2 <- rbind(resm2,res2)
		resm3 <- rbind(resm3,res3)
	}

	palette(gray(seq(0,.9,len=10)))
	op <- par(mfrow = c(1,3))
	plot(seq(-1,1,0.01),seq(0,1,0.005),col="blue",type='n')	
	for(y in 1:10){
		lines(seq(-1,1,0.01),resm[,y],col=palette()[y])
	}
	plot(seq(-1,1,0.01),seq(0,1,0.005),col="blue",type='n')	
	for(y in 1:10){
		lines(seq(-1,1,0.01),resm2[,y],col=palette()[y])
	}
	plot(seq(-1,1,0.01),seq(0,1.0,0.005),col="blue",type='n')	
	for(y in 1:10){
		lines(seq(-1,1,0.01),(resm3[,y]),col=palette()[y])
	}	
	
}

paper.plotHYPER <- function(){
	data(hyper)
	hyper <- fill.geno(hyper)
	cof <- MQMCofactorsEach(hyper,3)
	aa <- scanMQM(fill.geno(hyper),cof,step.size=2,step.max=150,step.min=0,windowsize=2)
	bb <- cim(hyper)
	cc <- scanone(hyper)
	plot(aa,bb,cc,chr=c(1,4,6),col=c("black","blue","red"),lwd=c(2,2,1),lty=c(1,6,1),main="Comparison QTL methodes: Dataset Hyper",ylab="LOD score")
	legend("topleft",c("MQM","CIM","MR"),col=c("black","blue","red"),lty=c(1,6,1),lwd=c(2,2,1))
}

paper.plotListeria <- function(){
	data(listeria)
	listeria <- fill.geno(listeria)
	cof <- MQMCofactorsEach(listeria,2)
	aa <- scanMQM(fill.geno(listeria),cof,step.size=2,step.max=100,step.min=0,windowsize=2)
	bb <- cim(listeria)
	cc <- scanone(listeria)
	plot(aa,bb,cc,chr=c(1,5,6,13),col=c("black","blue","red"),lwd=c(2,2,1),lty=c(1,6,1),main="Comparison QTL methodes: Dataset Listeria",ylab="LOD score")
	legend("left",c("MQM","CIM","MR"),col=c("black","blue","red"),lty=c(1,6,1),lwd=c(2,2,1))
}

paper.plotBristle3 <- function(){
	data(bristle3)
	bristle3 <- fill.geno(bristle3)
	cof <- MQMCofactorsEach(bristle3,3)
	aa <- scanMQMall(fill.geno(bristle3),cof,step.size=2,step.max=125,step.min=0,windowsize=2)
	bb <- cim(bristle3)
	cc <- scanone(bristle3)
	op <- par(mfrow = c(1,3))
	plot(aa,col=c("black"),lwd=c(2,2,1),lty=c(1,6,1),main="MQM: Dataset Bristle3",ylab="LOD score")
	plot(bb,col=c("blue"),lwd=c(2,2,1),lty=c(1,6,1),main="CIM: Dataset Bristle3",ylab="LOD score")
	plot(cc,col=c("red"),lwd=c(2,2,1),lty=c(1,6,1),main="MR: Dataset Bristle3",ylab="LOD score")
	aa <- scanMQMall(fill.geno(bristle3),cof,step.size=2,step.max=125,step.min=0,windowsize=2)
	plot.MQMnice(aa)
}

Lnormal <- function(res,vari){
  ret <- NULL
  ret <- exp(-0.5*(res/sqrt(vari))^2-log(sqrt(2*pi*vari)))
  cat("Var",vari,"Res",res,"->",ret,"\n")
  ret
}

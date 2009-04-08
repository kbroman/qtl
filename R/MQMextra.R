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


Lnormal <- function(res,vari){
  ret <- NULL
  ret <- exp(-0.5*(res/sqrt(vari))^2-log(sqrt(2*pi*vari)))
  cat("Var",vari,"Res",res,"->",ret,"\n")
  ret
}

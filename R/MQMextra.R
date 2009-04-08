extra <- function(a,b,x){
	.C("R_gammaln",as.double(3.0))
	lgamma(a)

	.C("R_betacf",as.double(3.0),as.double(3.0),as.double(0.9))

	.C("r_betai",a= as.double(3.0),b= as.double(3.0),x= as.double(0.9))
	pbeta(x,a,b)
	
	.C("R_Lnorm",as.double(3),as.double(8))
	
	resm <- NULL
	resm2 <- NULL
	for(x in seq(-1,1,0.01)){
		res <- NULL
		res2 <- NULL
		for(y in 1:2){
			aa <- .C("R_Lnorm",as.double(x),as.double(y))
			bb <- dnorm(x,0,y)
			res <- c(res,aa[2][[1]])
			res2 <- c(res2,bb)
		}
		resm <- rbind(resm,res)
		resm2 <- rbind(resm2,res2)
	}
	plot(seq(-1,1,0.01),resm[,1],col="blue")
	points(seq(-1,1,0.01),resm2[,1],col="green")
	points(seq(-1,1,0.01),resm[,2],col="lightblue")
	points(seq(-1,1,0.01),resm2[,2],col="lightgreen")

}
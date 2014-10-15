# scanonevar
# single-QTL genome scan for QTL affecting variance
# with code from Lars Ronnegard

vQTL <-
function( Rqtl.genoprob,pheno, vQTL.formula = ~add, vQTL.dformula = ~add,
         dglm.family = gaussian(), dglm.maxit = 25 , fixed.gamma.disp = FALSE)
{
    library(dglm)

    if (class(Rqtl.genoprob)[1]!="bc" & class(Rqtl.genoprob)[1]!="dh" & class(Rqtl.genoprob)[1]!="f2" )
        stop("This is not a backcross, double haploid, nor an F2")

    if (class(Rqtl.genoprob)[1]=="f2") {
        chrtype <- sapply(Rqtl.genoprob$geno, class)
        if (any(chrtype=="X")) Rqtl.genoprob <- subset( Rqtl.genoprob, chr=(chrtype != "X") )
    }

  formula.as.text<-paste(vQTL.formula,sep="")
  if (substr(formula.as.text[2],1,3)!="add")
      stop("Formula must be specified with add as first explanatory variable")
  formula.in<-as.formula(paste("pheno",formula.as.text[1],formula.as.text[2:length(formula.as.text)]))
  formula.as.text<-paste(vQTL.dformula,sep="")
  if (substr(formula.as.text[2],1,3)!="add")
      stop("Dformula must be specified with add as first explanatory variable")
  dformula.in<-as.formula(paste(formula.as.text[1],formula.as.text[2:length(formula.as.text)]))

###########
  ##BOX-COX PART
    library(MASS)
    formula.as.text<-paste(vQTL.formula,sep="")
	n.char <- nchar(formula.as.text[2])
	if (n.char>3) x.eff <- substr(formula.as.text[2],6, (n.char))
	if (n.char<4) x.eff <-"1"
	pheno.pos=pheno-(min(pheno)<0)*min(pheno)+0.1
	formula.bc <- as.formula(paste("pheno.pos",formula.as.text[1],x.eff))
	boxc <- boxcox(formula.bc,plotit=FALSE)
	#lambda <-boxc$x[sort(boxc$y,index.return=TRUE,decreasing=TRUE)$ix[1]]
	lambda <- boxc$x[which.max(boxc$y)]
    if (lambda<0.6 | lambda>1.6)
        print(paste("WARNING: Box-Cox transformation needed with lambda=", round(lambda,3)))
##########

  N=length(pheno) #No. of individuals
  n.chr=nchr(Rqtl.genoprob) #No. of chromosomes
  n.plotrows=1+(n.chr>5)+(n.chr>10)+(n.chr>15)+(n.chr>20)
  par(mfrow=c(n.plotrows,ceiling(n.chr/n.plotrows)))
  chr.names <- names(Rqtl.genoprob$geno)
  out <- matrix(0, n.chr, 4)
  colnames(out) = c("max(logP.mean)", "position", "max(logP.dispersion)", "position")
  rownames(out) = chr.names
  scan.logPm <- scan.logPd <- chr.names.out <- NULL
  for (j in 1:n.chr){
    max.logP.m <- max.logP.d  <- 0
    max.pos.m <- max.pos.d <- NA
    if (class(Rqtl.genoprob)[1]=="f2") {
      g11<-Rqtl.genoprob$geno[[j]]$prob[1:N,,1]
      g12<-Rqtl.genoprob$geno[[j]]$prob[1:N,,2]
      g13<-Rqtl.genoprob$geno[[j]]$prob[1:N,,3]
      a1<-g11+0.5*g12
    }
    else {
      a1<-Rqtl.genoprob$geno[[j]]$prob[1:N,,1]
    }
    xlab.text=paste("Chr",chr.names[j])
    n.loci=dim(a1)[2]
	 logP.m<-logP.d<-numeric(n.loci)
	 for (i in 1:n.loci) {
		add<-as.numeric(a1[,i])
		d.fit<-try(dglm(formula=formula.in,dformula=dformula.in,
                        family = dglm.family, control=dglm.control(maxit=dglm.maxit)),silent=TRUE)
    if ((class(d.fit)!="dglm")[1]) {
			print(paste("ERROR: dglm did not converge on chromosome",chr.names[j]," position",i))
			print("CHECK MARKER INFORMATION")
			}
		else {
			p.mean<-summary(d.fit)$coef[2,4]
			if (!fixed.gamma.disp) p.disp<- summary(d.fit$dispersion)$coef[2,4]
			if (fixed.gamma.disp) p.disp<-(summary(d.fit)$dispersion.summary$coeff[2,4])
			logP.m[i]<- -log10(p.mean)
			logP.d[i]<- -log10(p.disp)
			if (-log10(p.mean) > max.logP.m) {
				max.logP.m = -log10(p.mean)
				max.pos.m = i
			}
			if (-log10(p.disp) > max.logP.d) {
				max.logP.d = -log10(p.disp)
				max.pos.d = i
			}
			if (d.fit$iter==dglm.maxit) {
				logP.d[i]=0
				print(paste("ERROR: dglm did not converge on chromosome",j," position",i))
				print("CHECK MARKER INFORMATION")
			}
		}
	 }
 	if (max(logP.m)<10 & max(logP.d)<10) {
		plot(logP.m,type="l",ylim=c(0,10),xlab=xlab.text,ylab="logP")
  		points(logP.d,type="l",col="red")
	}
	else {
		plot(logP.m,type="l",ylim=c(0,max(c(logP.m,logP.d))),xlab=xlab.text,ylab="logP")
		abline(10,0,lty=2)
  		points(logP.d,type="l",col="red")
	}
	out[j,1:4] <- c(max.logP.m,max.pos.m, max.logP.d,max.pos.d)
	#chr.names.out <- c(chr.names.out, rep(n.loci,chr.names[j]))
	scan.logPm <- c(scan.logPm, logP.m)
	scan.logPd <- c(scan.logPd, logP.d)
  }
  print(out)

  list(QTL.table = out, scan.QTL = scan.logPm, scan.vQTL = scan.logPd, Rqtl.obj = Rqtl.genoprob)
}


# scanonevar
# single-QTL genome scan for QTL affecting variance
# with code from Lars Ronnegard

scanonevar <-
function(cross, pheno.col=1, mean_formula = ~add, var_formula = ~add,
         family = gaussian(), maxit = 25 , fixed_gamma_disp = FALSE,
         quiet=TRUE)
{
    # check input
    crosstype <- class(cross)[1]
    if(!(crosstype %in% c("bc", "dh", "f2", "haploid", "risib", "riself")))
       stop('scanonevar not implemented for cross type "', crosstype, '"')
    chrtype <- sapply(cross$geno, class)
    if(any(chrtype=="X")) {
        warning("Analysis of X chromosome not implemented for scanonevar; omitted.")
        cross <- subset(cross, chr=(chrtype != "X"))
    }

    # grab phenotype
    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }
    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(any(is.na(num))) {
            if(sum(is.na(num)) > 1)
                stop("Couldn't identify phenotypes ", paste(paste("\"", pheno.col[is.na(num)], "\"", sep=""),
                                                            collapse=" "))
            else
                stop("Couldn't identify phenotype \"", pheno.col[is.na(num)], "\"")
        }
        pheno.col <- num
    }
    if(any(pheno.col < 1 | pheno.col > nphe(cross)))
        stop("pheno.col values should be between 1 and the no. phenotypes")
    pheno <- cross$pheno[,pheno.col]
    if(is.matrix(pheno) && ncol(pheno) > 1) {
        pheno <- pheno[,1]
        warning('scanonevar requires a single phenotype; all but "', phenames(cross)[pheno.col[1]], '" omitted.')
    }

    # mean and variance formulas
    formula.as.text<-paste(mean_formula,sep="")
    if(substr(formula.as.text[2],1,3)!="add")
        stop("Formula must be specified with add as first explanatory variable")
    formula.in <- as.formula(paste("pheno", formula.as.text[1], formula.as.text[2:length(formula.as.text)]))
    formula.as.text<-paste(var_formula, sep="")
    if(substr(formula.as.text[2],1,3)!="add")
        stop("Dformula must be specified with add as first explanatory variable")
    dformula.in <- as.formula(paste(formula.as.text[1], formula.as.text[2:length(formula.as.text)]))

    ###########
    ##BOX-COX PART
    formula.as.text<-paste(var_formula,sep="")
	n.char <- nchar(formula.as.text[2])
	if (n.char>3) x.eff <- substr(formula.as.text[2],6, (n.char))
	if (n.char<4) x.eff <-"1"
	pheno.pos <- pheno-(min(pheno)<0)*min(pheno)+0.1
	formula.bc <- as.formula(paste("pheno.pos",formula.as.text[1],x.eff))
	boxc <- MASS::boxcox(formula.bc,plotit=FALSE)
	#lambda <-boxc$x[sort(boxc$y,index.return=TRUE,decreasing=TRUE)$ix[1]]
	lambda <- boxc$x[which.max(boxc$y)]
    if (lambda<0.6 | lambda>1.6)
        warning("Box-Cox transformation needed with lambda = ", round(lambda,3))
    ##########

    N <- length(pheno) #No. of individuals
    n.chr <- nchr(cross) #No. of chromosomes
    chr.names <- chrnames(cross)

    # need to run calc.genoprob?
    if(!("prob" %in% names(cross$geno[[1]]))) {
        warning("First running calc.genoprob")
        cross <- calc.genoprob(cross)
    }

    scan.logPm <- scan.logPd <- chr.names.out <- NULL

    result <- NULL
    for(j in seq(along=cross$geno)) {
        if(!quiet) message("Chr ", chr.names[j])

        if (crosstype=="f2") {
            g11 <- cross$geno[[j]]$prob[,,1]
            g12 <- cross$geno[[j]]$prob[,,2]
            a1  <- g11 + g12/2
        }
        else {
            a1 <- cross$geno[[j]]$prob[,,1]
        }

        n.loci <- dim(a1)[2]

        logP.m <- logP.d <- numeric(n.loci)

        for(i in 1:n.loci) {
            add <- as.numeric(a1[,i])
            d.fit<-try(dglm::dglm(formula=formula.in,dformula=dformula.in,
                                  family = family, control=dglm::dglm.control(maxit=maxit)),silent=TRUE)

            if ((class(d.fit)!="dglm")[1]) {
                warning("dglm did not converge on chr ", chr.names[j], " position ", i)
			}
            else {
                p.mean <- summary(d.fit)$coef[2,4]
                if (!fixed_gamma_disp)
                    p.disp<- summary(d.fit$dispersion)$coef[2,4]
                if (fixed_gamma_disp)
                    p.disp<-(summary(d.fit)$dispersion.summary$coeff[2,4])
                logP.m[i]<- -log10(p.mean)
                logP.d[i]<- -log10(p.disp)
                if (d.fit$iter==maxit) {
                    logP.d[i]=0
                    warning("dglm did not converge on chr", chr.names[j], " position ", i)
                }
            }
        }

        map <- attr(cross$geno[[j]]$prob,"map")
        w <- names(map)
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",chr.names[j],".",w[o],sep="")
        thischr <- data.frame(chr=rep(chr.names[j], length(w)),
                              pos=map,
                              neglogP_mean=logP.m,
                              neglogP_disp=logP.d, stringsAsFactors=FALSE)
        rownames(thischr) <- w

        if(is.null(result)) result <- thischr
        else result <- rbind(result, thischr)
    }

    class(result) <- c("scanone", "data.frame")
    attr(result, "method") <- "scanonevar"

    result
}


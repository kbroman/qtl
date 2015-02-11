# scanonevar.meanperm
# single-QTL genome scan for QTL affecting variance
# with code from Lars Ronnegard


scanonevar.meanperm <-
  function(cross, pheno.col=1, mean_covar = NULL, var_covar = NULL,
           maxit = 25 , tol=1e-6, n.mean.perm = 2, seed = 27517, quiet=TRUE)
  {

    set.seed(seed)

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
    #     if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
    #       cross$pheno <- cbind(pheno.col, cross$pheno)
    #       pheno.col <- 1
    #     }
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

    N <- length(pheno) # No. individuals
    n.chr <- nchr(cross) #No. chromosomes
    chr.names <- chrnames(cross)

    # need to run calc.genoprob?
    if(!("prob" %in% names(cross$geno[[1]]))) {
      warning("First running calc.genoprob")
      cross <- calc.genoprob(cross)
    }

    scan.logPm <- scan.logPd <- chr.names.out <- NULL

    # set up data and formulas
    X <- cbind(pheno=pheno,
               mean.add=rep(0, length(pheno)),
               var.add=rep(0, length(pheno)))
    mean_formula <- "pheno ~ mean.add"
    var_formula <- "~ var.add"

    if(!is.null(mean_covar)) {
      ncolX <- ncol(X)
      X <- cbind(X, mean_covar)
      meancovarnames <- paste0("meancov", 1:(ncol(X)-ncolX))
      colnames(X)[-(1:ncolX)] <- meancovarnames
      mean_formula <- paste(mean_formula, "+",
                            paste(meancovarnames, collapse="+"))
    }
    if(!is.null(var_covar)) {
      ncolX <- ncol(X)
      X <- cbind(X, var_covar)
      varcovarnames <- paste0("varcov", 1:(ncol(X)-ncolX))
      colnames(X)[-(1:ncolX)] <- varcovarnames
      var_formula <- paste(var_formula, "+",
                           paste(varcovarnames, collapse="+"))
    }
    X <- as.data.frame(X)
    mean_formula <- as.formula(mean_formula)
    var_formula <- as.formula(var_formula)

    max.mean.neglog.ps <- rep(NA, n.mean.perm)

    for(perm.num in 1:n.mean.perm) {

      result <- NULL
      for(j in seq(along=cross$geno)) { # loop over chromosomes
        if(!quiet) message(" - Chr ", chr.names[j])

        if (crosstype=="f2") {
          g11 <- cross$geno[[j]]$prob[,,1]
          g12 <- cross$geno[[j]]$prob[,,2]
          g13 <- cross$geno[[j]]$prob[,,3]
          a1  <- g11 + g12/2
          d1 <-  g12 - (g11+g13)/2
        }
        else {
          a1 <- cross$geno[[j]]$prob[,,1]
        }

        n.loci <- dim(a1)[2]

        logP.m <- logP.d <- numeric(n.loci)

        for(i in 1:n.loci) { # loop over positions within chromosome
          # fill in genotype probs for this locus
          X$mean.add <- sample(a1[,i])
          X$var.add <- a1[,i]

          d.fit <- DGLM_norm(m.form=mean_formula, d.form=var_formula, indata=X,
                             maxiter=maxit, conv=tol)

          p.mean <- summary(d.fit$mean)$coef[2,4]
          p.disp<- summary(d.fit$disp)$coef[2,4]
          if (d.fit$iter < maxit) {
            logP.m[i]<- -log10(p.mean)
            logP.d[i]<- -log10(p.disp)
          }
          else {
            logP.m[i]<- -log10(p.mean)
            logP.d[i]<- 0
            warning("dglm did not converge on chr", chr.names[j], " position ", i)
          }
        }

        max.mean.neglog.ps[perm.num] <- max(max.mean.neglog.ps[perm.num], logP.m, na.rm = TRUE)
      }

    }

    return(max.mean.neglog.ps)
  }

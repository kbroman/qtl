
######################################################################
# count terms in a model, for use by plod, Quoc editted, intermediate function to get onX
######################################################################
countqtltermsQuoc <-
  function(formula, qtl, ignore.covar=TRUE) # qtl is used to extract a logical vector: on "X" chr == 1
  {
    if(is.character(formula)) formula <- as.formula(formula)
    factors <- attr(terms(formula), "factors")[-1,,drop=FALSE]
    if(any(factors > 1))  {
      warning("some formula terms > 1; may be a problem with the formula:\n    ", deparseQTLformula(formula))
      factors[factors > 1] <- 1
    }
    nterm <- apply(factors, 2, sum)
    
    if(any(nterm>2)) 
      stop("Can't deal with higher-order interactions\n")
    
    # need to check for QTL x covariate interactions in here!
    
    if(ignore.covar) {
      cn <- colnames(factors)
      wh <- c(grep("^[Qq][0-9]+$", cn),
              grep("^[Qq][0-9]+:[Qq][0-9]+$", cn))
      rn <- rownames(factors)
      wh2 <- c(grep("^[Qq][0-9]+$", rn),
               grep("^[Qq][0-9]+:[Qq][0-9]+$", rn))
      factors <- factors[wh2,wh, drop=FALSE]
    }
    nterm <- apply(factors, 2, sum)
    
    nmain <- sum(nterm==1)
    
    qtltype <- qtl$chrtype
    names(qtltype) <- qtl$name
    factors.main <- factors[,nterm==1, drop=FALSE] # Get first order factors
    nmainA <- sum(qtltype[colnames(factors.main)]=="A")
    nmainX <- sum(qtltype[colnames(factors.main)]=="X")

    
    if(all(nterm==1))
      return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=0, intAX=0, intXX=0, inttot=0))
    
    n.int <- sum(nterm==2)
    
    if(n.int <=0) # 0 interactions, so no need to figure them out
      return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=n.int, intAX=0, intXX=0, inttot=n.int))
    
    factors <- factors[, nterm==2, drop=FALSE] # Get second order factors
    
    wh <- apply(factors, 2, function(a) which(a==1)) # Get the location of interaction, check for X here.
        
    if(n.int ==1){ # 1 interaction, check if it is AA AX or XX
      int.type <- sum(qtltype[rownames(factors)[wh[,1]]]=="X") #0 is AA, 1 is AX, 2 is XX
      if (int.type == 0) return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=1, intAX=0, intXX=0, inttot=n.int))
      else 
        if (int.type == 1) return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=0, intAX=1, intXX=0, inttot=n.int))
        else return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=0, intAX=0, intXX=1, inttot=n.int))
    }
    
    u <- sort(unique(as.numeric(wh))) # Unique nodes with interactions in increasing order
    grp <- rep(NA, length(u)) # Group member, length = number of node
    names(grp) <- u
    
    ngrp <- 0 # Number of group
    nintAA <- NULL # Number of AA interaction of specific group
    nintAX <- NULL # Number of AX interaction of specific group
    nintXX <- NULL # Number of XX interaction of specific group
    
    
    for(i in 1:ncol(wh)) {
      thegrp <- grp[as.character(wh[,i])]
      int.type <- sum(qtltype[rownames(factors)[wh[,i]]]=="X") #0 is AA, 1 is AX, 2 is XX
      if(all(!is.na(thegrp))) { # Merge 2 groups
        nintAA[as.character(thegrp[1])] <-
          sum(nintAA[unique(as.character(thegrp))])
        nintAX[as.character(thegrp[1])] <-
          sum(nintAX[unique(as.character(thegrp))])
        nintXX[as.character(thegrp[1])] <-
          sum(nintXX[unique(as.character(thegrp))])
        if (int.type==0) nintAA[as.character(thegrp[1])] <- nintAA[as.character(thegrp[1])] + 1
        if (int.type==1) nintAX[as.character(thegrp[1])] <- nintAX[as.character(thegrp[1])] + 1
        if (int.type==2) nintXX[as.character(thegrp[1])] <- nintXX[as.character(thegrp[1])] + 1
        grp[grp==thegrp[1] | grp==thegrp[2]] <- thegrp[1]
      } # two connected group become one group and has the name of the first group, number of interaction is updated
      else if(any(!is.na(thegrp))) { # add 1 more interaction to current group
        grp[as.character(wh[,i])] <- thegrp[!is.na(thegrp)]
        if (int.type==0) nintAA[as.character(thegrp[!is.na(thegrp)])] <- 
          nintAA[as.character(thegrp[!is.na(thegrp)])] + 1
        if (int.type==1) nintAX[as.character(thegrp[!is.na(thegrp)])] <- 
          nintAX[as.character(thegrp[!is.na(thegrp)])] + 1
        if (int.type==2) nintXX[as.character(thegrp[!is.na(thegrp)])] <- 
          nintXX[as.character(thegrp[!is.na(thegrp)])] + 1
      }
      else { # introduce new group
        ngrp <- ngrp+1
        grp[as.character(wh[,i])] <- ngrp
        # Initialize nint
        nintAA[as.character(ngrp)] <- 0
        nintAX[as.character(ngrp)] <- 0
        nintXX[as.character(ngrp)] <- 0
        if (int.type==0) nintAA[as.character(ngrp)] <- 1
        if (int.type==1) nintAX[as.character(ngrp)] <- 1
        if (int.type==2) nintXX[as.character(ngrp)] <- 1
      }
    }
    
    nintAA <- nintAA[as.character(unique(grp))]
    nintAX <- nintAX[as.character(unique(grp))]
    nintXX <- nintXX[as.character(unique(grp))]
    nL <- sum(nintAA>0)
    nH <- sum(nintAA)-nL
    c(mainA=nmainA, mainX=nmainX, intH=nH, intL=nL, intAX=sum(nintAX), intXX=sum(nintXX), inttot=nH+nL+sum(nintAX)+sum(nintXX))
  }


######################################################################
# calculate penalties for pLOD using scantwo permutation results,  Quoc editted
######################################################################
calcQuoc.penalties <-
  function(perms, alpha=0.05, lodcolumn)
  {
    if(missing(perms) || !("scantwoperm" %in% class(perms)))
      stop("You must include permutation results from scantwo.")
    
    if(missing(lodcolumn)) {
      if(is.matrix(perms[[1]]) && ncol(perms[[1]]) > 1)
        lodcolumn <- 1:ncol(perms[[1]])
      else lodcolumn <- 1
    }
    
    if(length(lodcolumn)>1) {
      result <- NULL
      for(i in seq(along=lodcolumn)) {
        temp <- calc.penalties(perms, alpha, lodcolumn[i])
        result <- rbind(result, temp)
      }
      dimnames(result) <- list(colnames(perms[[1]])[lodcolumn], names(temp))
      return(result)
    }
    
    if(is.matrix(perms[[1]]) && ncol(perms[[1]]) >1) {
      if(lodcolumn < 1 || lodcolumn > ncol(perms[[1]]))
        stop("lodcolumn misspecified")
      for(i in seq(along=perms))
        perms[[i]] <- perms[[i]][,lodcolumn,drop=FALSE]
    }
    
    qu <- summary(perms, alpha=alpha)
    if(!("one" %in% names(qu)))
      stop("You need to re-run scantwo permutations with R/qtl version >= 1.09.")
    
    qu <- summary(perms, alpha=alpha)
    
    if(length(alpha)>1) {
      penalties <- cbind(qu[["one"]], qu[["trueint"]], qu[["fv1"]]-qu[["one"]])
      colnames(penalties) <- c("main","heavy", "light")
    }
    else {
      penalties <- c(qu[["one"]], qu[["trueint"]], qu[["fv1"]]-qu[["one"]])
      names(penalties) <- c("main","heavy", "light")
    }
    penalties
  }
######################################################################
stepwiseqtlQuoc <-
  function(cross, chr, pheno.col=1, qtl, formula, max.qtl=10, covar=NULL,
           method=c("imp", "hk"), model=c("normal", "binary"), incl.markers=TRUE, refine.locations=TRUE,
           additive.only=FALSE, scan.pairs=FALSE, penalties,
           keeplodprofile=FALSE, keeptrace=FALSE, verbose=TRUE,
           tol=1e-4, maxit=1000)
  {
    if(!("cross" %in% class(cross)))
      stop("Input should have class \"cross\".")
    
    if(!missing(chr)) cross <- subset(cross, chr)
    
    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
      cross$pheno <- cbind(pheno.col, cross$pheno)
      pheno.col <- 1
    }
    
    if(!missing(qtl)) { # start f.s. at somewhere other than the null
      if( !("qtl" %in% class(qtl)) )
        stop("The qtl argument must be an object of class \"qtl\".")
      
      # check that chromosomes were retained, otherwise give error
      m <- is.na(match(qtl$chr, names(cross$geno)))
      if(any(m)) {
        wh <- qtl$chr[m]
        if(length(wh) > 1)
          stop("Chromosomes ", paste(wh, collapse=", "), " (in QTL object) not in cross object.")
        else
          stop("Chromosome ", wh, " (in QTL object) not in cross object.")
      }
      if(missing(formula)) { # create a formula with all covariates and all QTL add've
        if(!is.null(covar)) 
          formula <- paste("y ~ ", paste(names(covar), collapse="+"), "+")
        else
          formula <- "y ~ "
        formula <- paste(formula, paste(paste("Q", 1:length(qtl$chr), sep=""), collapse="+"))
      }
      else {
        temp <- checkStepwiseqtlStart(qtl, formula, covar)
        qtl <- temp$qtl
        formula <- temp$formula
      }
      startatnull <- FALSE
    }
    else {
      if(!missing(formula)) 
        warning("formula ignored if qtl is not provided.")
      startatnull <- TRUE
    }
    
    # revise names in qtl object
    if(!startatnull) 
      qtl$name <- qtl$altname
    
    # check that we have the right stuff for the selected method
    method <- match.arg(method)
    model <- match.arg(model)
    if(method=="imp") {
      if(!("draws" %in% names(cross$geno[[1]]))) {
        if("prob" %in% names(cross$geno[[1]])) {
          warning("The cross doesn't contain imputations; using method=\"hk\".")
          method <- "hk"
        }
        else
          stop("You need to first run sim.geno.")
      }
    }
    else {
      if(!("prob" %in% names(cross$geno[[1]]))) {
        if("draws" %in% names(cross$geno[[1]])) {
          warning("The cross doesn't contain QTL genotype probabilities; using method=\"imp\".")
          method <- "imp"
        }
        else
          stop("You need to first run calc.genoprob.")
      }
    }
    if(method=="imp") qtlmethod <- "draws"
    else qtlmethod <- "prob"
    
    if(!missing(qtl) && qtl$n.ind != nind(cross)) {
      warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
      if(method=="imp")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
      else
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="prob")
    }
    
    if(!missing(qtl) && method=="imp" && dim(qtl$geno)[3] != dim(cross$geno[[1]]$draws)[3])  {
      warning("No. imputations in qtl object doesn't match that in the input cross; re-creating qtl object.")
      qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
    }    
    
    # check that qtl object matches the method
    if(!startatnull) {
      if(method=="imp" && !("geno" %in% names(qtl)))
        stop("The qtl object doesn't contain imputations; re-run makeqtl with what=\"draws\".")
      else if(method=="hk" && !("prob" %in% names(qtl))) 
        stop("The qtl object doesn't contain QTL genotype probabilities; re-run makeqtl with what=\"prob\".")
    }
    
    # check phenotypes and covariates; drop ind'ls with missing values
    if(length(pheno.col) > 1) {
      pheno.col <- pheno.col[1]
      warning("stepwiseqtl can take just one phenotype; only the first will be used")
    }
    if(is.character(pheno.col)) {
      num <- find.pheno(cross, pheno.col)
      if(is.na(num)) 
        stop("Couldn't identify phenotype \"", pheno.col, "\"")
      pheno.col <- num
    }
    if(any(pheno.col < 1 | pheno.col > nphe(cross)))
      stop("pheno.col values should be between 1 and the no. phenotypes")
    pheno <- cross$pheno[,pheno.col]
    if(!is.null(covar)) phcovar <- cbind(pheno, covar)
    else phcovar <- as.data.frame(pheno, stringsAsFactors=TRUE)
    hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
    if(all(hasmissing))
      stop("All individuals are missing phenotypes or covariates.")
    if(any(hasmissing)) {
      pheno <- pheno[!hasmissing]
      cross <- subset(cross, ind=!hasmissing)
      if(!is.null(covar)) covar <- covar[!hasmissing,,drop=FALSE]
      
      if(!startatnull) {
        if(method=="imp") 
          qtl$geno <- qtl$geno[!hasmissing,,,drop=FALSE]
        else {
          for(i in seq(along=qtl$prob))
            qtl$prob[[i]] <- qtl$prob[[i]][!hasmissing,,drop=FALSE]
        }
        qtl$n.ind <- sum(!hasmissing)
      }
    }
    
    if(max.qtl < 1)
      stop("Need max.qtl > 0 if we are to scan for qtl")
    
    if(is.null(covar)) {
      lod0 <- 0
      if(startatnull)
        firstformula <- y~Q1
      else firstformula <- formula
    }
    else {
      lod0 <- length(pheno)/2 * log10(sum((pheno-mean(pheno))^2) / sum(lm(pheno ~ as.matrix(covar))$resid^2))
      if(startatnull)
        firstformula <- as.formula(paste("y~", paste(names(covar), collapse="+"), "+", "Q1"))
      else firstformula <- formula
    }
    
    
    # penalties, Quoc: need change later, comment out for now
    cross.type <- class(cross)[1]
#     if(missing(penalties)) {
#       if(cross.type=="f2") {
#         penalties <-  c(3.52, 4.28, 2.69)
#       }
#       else if(cross.type=="bc") {
#         penalties <-  c(2.69, 2.62, 1.19)
#       }
#       else 
#         stop("No default penalties available for cross type ", cross.type)
#     }
#     else if(length(penalties) != 3) {
#       if(length(penalties)==1) {
#         if(additive.only)
#           penalties <- c(penalties,Inf,Inf)
#         else
#           stop("You must include a penalty for interaction terms.")
#       }
#       else {
#         if(length(penalties)==2) 
#           penalties <- penalties[c(1,2,2)]
#         else {
#           warning("penalties should have length 3")
#           penalties <- penalties[1:3]
#         }
#       }
#     }
#     
    if(verbose > 2) verbose.scan <- TRUE
    else verbose.scan <- FALSE
    
    curbest <- NULL
    curbestplod <- 0
    
    # initial scan : either 1d or 2d
    if(verbose) cat(" -Initial scan\n")
    if(startatnull) { # Quoc: changed
      if(additive.only || max.qtl == 1 || !scan.pairs) {
        out <- scanone(cross, pheno.col=pheno.col, method=method, model="normal",
                       addcovar=covar)
        lod <- max(out[,3], na.rm=TRUE)
        
        
        wh <- which(!is.na(out[,3]) & out[,3]==lod)
        if(length(wh) > 1) wh <- sample(wh, 1)
        
        qtl <- makeqtl(cross, as.character(out[wh,1]), out[wh,2], "Q1",
                       what=qtlmethod)
        formula <- firstformula
        n.qtl <- 1
        curplod <- calc.plodQuoc(lod, formula=firstformula, qtl=qtl, penalties=penalties) # update plod after choosing formula and qtl
      }
      else {
        out <- scantwo(cross, pheno.col=pheno.col, method=method, model="normal",
                       incl.markers=incl.markers, addcovar=covar, verbose=verbose.scan)
        lod <- out$lod
        #pLOD 1 term
        lod1 <- max(diag(lod), na.rm=TRUE)
        wh <- which(!is.na(diag(lod)) & diag(lod) == lod1)
        if(length(wh) > 1) wh <- sample(wh, 1)
        m <- out$map[wh,]
        qtl1 <- makeqtl(cross, as.character(m[1,1]), m[1,2], "Q1", what=qtlmethod)
        formula1 <- firstformula
        plod1 <- calc.plodQuoc(lod1, qtl=qtl1, formula=formula1, penalties=penalties)
        #pLOD 2 additive terms        
        loda <- max(lod[upper.tri(lod)], na.rm=TRUE)
        temp <- max(out, what="add")
        if(nrow(temp) > 1)
          temp <- temp[sample(1:nrow(temp),1),]
        qtla <- makeqtl(cross, c(as.character(temp[1,1]), as.character(temp[1,2])),
                       c(temp[1,3], temp[1,4]), c("Q1","Q2"), what=qtlmethod)
        formulaa <- as.formula(paste(deparseQTLformula(firstformula), "+Q2", sep=""))
        ploda <- calc.plod(loda, qtl=qtla, formula=formulaa,
                           penalties=penalties)
        #pLOD 2 terms with interaction
        lodf <- max(lod[lower.tri(lod)], na.rm=TRUE)
        temp <- max(out, what="full")
        if(nrow(temp) > 1)
          temp <- temp[sample(1:nrow(temp),1),]
        qtlf <- makeqtl(cross, c(as.character(temp[1,1]), as.character(temp[1,2])),
                       c(temp[1,3], temp[1,4]), c("Q1","Q2"), what=qtlmethod)
        formulaf <- as.formula(paste(deparseQTLformula(firstformula), "+Q2+Q1:Q2", sep=""))
        plodf <- calc.plod(lodf, qtl=qtlf, formula=formulaf,
                           penalties=penalties)
        
        if(plod1 > ploda && plod1 > plodf) {
          qtl <- qtl1
          formula <- formula1
          curplod <- plod1
          lod <- lod1
          n.qtl <- 1
        }
        else if(ploda > plodf) {
          qtl <- qtla
          formula <- formulaa
          curplod <- ploda
          lod <- loda
          n.qtl <- 2
        }
        else {
          qtl <- qtlf
          formula <- formulaf
          curplod <- plodf
          lod <- lodf
          n.qtl <- 2
        }
      }
    } # start at null, Quoc: changed
    else {
      if(verbose) cat(" ---Starting at a model with", length(qtl$chr), "QTL\n")
      if(refine.locations) {
        if(verbose) cat(" ---Refining positions\n")
        rqtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl,
                          covar=covar, formula=formula, method=method,
                          verbose=verbose.scan, incl.markers=incl.markers,
                          keeplodprofile=FALSE)
        if(any(rqtl$pos != qtl$pos)) { # updated positions
          if(verbose) cat(" ---  Moved a bit\n")
        }
        qtl <- rqtl
      }
      lod <- fitqtl(cross, pheno.col, qtl, covar=covar, formula=formula,
                    method=method, model=model, dropone=FALSE, get.ests=FALSE,
                    run.checks=FALSE, tol=tol, maxit=maxit)$result.full[1,4] - lod0
      curplod <- calc.plodQuoc(lod, formula=formula, qtl=qtl,
                           penalties=penalties) # Quoc changed
      attr(qtl, "pLOD") <- curplod
      n.qtl <- length(qtl$chr)
    }
    
    attr(qtl, "formula") <- deparseQTLformula(formula)
    attr(qtl, "pLOD") <- curplod
    
    if(curplod > 0) {
      curbest <- qtl
      curbestplod <- curplod
      
      if(verbose) 
        cat("** new best ** (pLOD increased by ", round(curplod, 4), ")\n", sep="")
    }
    
    if(keeptrace) {
      temp <- list(chr=qtl$chr, pos=qtl$pos)
      attr(temp, "formula") <- deparseQTLformula(formula)
      attr(temp, "pLOD") <- curplod
      class(temp) <- c("compactqtl", "list")
      thetrace <- list("0"=temp)
    }
    
    if(verbose)
      cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
          deparseQTLformula(formula), "\n")
    if(verbose > 1) 
      cat("         qtl:", paste(qtl$chr, round(qtl$pos,1), sep="@"), "\n")
    
    # start stepwise search
    i <- 0
    while(n.qtl < max.qtl) {
      i <- i+1
      
      if(verbose) {
        cat(" -Step", i, "\n")
        cat(" ---Scanning for additive qtl\n")
      }
      
      out <- addqtl(cross, pheno.col=pheno.col, qtl=qtl, covar=covar,
                    formula=formula, method=method, incl.markers=incl.markers,
                    verbose=verbose.scan)
      
      curlod <- max(out[,3], na.rm=TRUE) 
      wh <- which(!is.na(out[,3]) & out[,3]==curlod)
      if(length(wh) > 1) wh <- sample(wh,1)
      curqtl <- addtoqtl(cross, qtl, as.character(out[wh,1]), out[wh,2],
                         paste("Q", n.qtl+1, sep=""))
      curformula <- as.formula(paste(deparseQTLformula(formula), "+Q", n.qtl+1, sep=""))
      curlod <- curlod + lod # Quoc: Why? Because LOD scores returned are relative to the base model  
      curplod <- calc.plodQuoc(curlod, formula=curformula, qtl=curqtl,
                           penalties=penalties)
      if(verbose) cat("        plod =", curplod, "\n")
      
      curnqtl <- n.qtl+1
      
      if(!additive.only) {
        for(j in 1:n.qtl) {
          
          if(verbose)
            cat(" ---Scanning for QTL interacting with Q", j, "\n", sep="")
          
          thisformula <- as.formula(paste(deparseQTLformula(formula), "+Q", n.qtl+1,
                                          "+Q", j, ":Q", n.qtl+1, sep=""))
          out <- addqtl(cross, pheno.col=pheno.col, qtl=qtl, covar=covar,
                        formula=thisformula, method=method, incl.markers=incl.markers,
                        verbose=verbose.scan)
          thislod <- max(out[,3], na.rm=TRUE)
          
          wh <- which(!is.na(out[,3]) & out[,3]==thislod)
          if(length(wh) > 1) wh <- sample(wh,1)
          thisqtl <- addtoqtl(cross, qtl, as.character(out[wh,1]), out[wh,2],
                              paste("Q", n.qtl+1, sep=""))
          
          thislod <- thislod + lod
          thisplod <- calc.plodQuoc(thislod, formula=thisformula, qtl=thisqtl,
                                penalties=penalties)
          if(verbose) cat("        plod =", thisplod, "\n")
          
          if(thisplod > curplod) {
            curformula <- thisformula
            curplod <- thisplod
            curlod <- thislod
            curqtl <- thisqtl
            
            curnqtl <- n.qtl+1
          }
        }
        
        if(n.qtl > 1) {
          if(verbose)
            cat(" ---Look for additional interactions\n")
          temp <- addint(cross, pheno.col, qtl, covar=covar, formula=formula,
                         method=method, qtl.only=TRUE, verbose=verbose.scan)
          if(!is.null(temp)) {
            thislod <- max(temp[,3], na.rm=TRUE)
            wh <- which(!is.na(temp[,3]) & temp[,3] == thislod)
            if(length(wh) > 1) wh <- sample(wh, 1)
            thisformula <- as.formula(paste(deparseQTLformula(formula), "+", rownames(temp)[wh]))
            thislod <- thislod + lod
            thisplod <- calc.plodQuoc(thislod, formula=thisformula, qtl=qtl,
                                  penalties=penalties)
            if(verbose) cat("        plod =", thisplod, "\n")
            if(thisplod > curplod) {
              curformula <- thisformula
              curplod <- thisplod
              curlod <- thislod
              curqtl <- qtl
              curnqtl <- n.qtl
            }
          }
        }
        
        if(scan.pairs) { #Quoc: changed, result always have even number of QTLs?
          if(verbose)
            cat(" ---Scan for an additional pair\n")
          out <- addpair(cross, pheno.col=pheno.col, qtl=qtl, covar=covar,
                         formula=formula, method=method, incl.markers=incl.markers,
                         verbose=verbose.scan)
          thelod <- out$lod
          
          loda <- max(thelod[upper.tri(thelod)], na.rm=TRUE)
          temp <- max(out, what="add")
          if(nrow(temp) > 1)
            temp <- temp[sample(1:nrow(temp),1),]
          qtlspa <- addtoqtl(cross, qtl, c(as.character(temp[1,1]), as.character(temp[1,2])),
                             c(temp[1,3], temp[1,4]), paste("Q", n.qtl+1:2, sep=""))
          formulaspa <- as.formula(paste(deparseQTLformula(formula), "+Q", n.qtl+1, "+Q",
                                         n.qtl+2, sep=""))
          ploda <- calc.plodQuoc(loda+lod, formula=formulaspa, qtl=qtlspa,
                             penalties=penalties) 
          lodf <- max(thelod[lower.tri(thelod)], na.rm=TRUE)
          temp <- max(out, what="full")
          if(nrow(temp) > 1)
            temp <- temp[sample(1:nrow(temp),1),]
          
          qtlspf <- addtoqtl(cross, qtl, c(as.character(temp[1,1]), as.character(temp[1,2])),
                             c(temp[1,3], temp[1,4]), paste("Q", n.qtl+1:2, sep=""))
          formulaspf <- as.formula(paste(deparseQTLformula(formula), "+Q", n.qtl+1, "+Q",
                                         n.qtl+2, "+Q", n.qtl+1, ":Q", n.qtl+2, 
                                         sep=""))
          plodf <- calc.plodQuoc(lodf+lod, formula=formulaspf, qtl=qtlspf,
                             penalties=penalties)
          
          if(verbose) {
            cat("        ploda =", ploda, "\n")
            cat("        plodf =", plodf, "\n")
          }
          
          if(ploda > curplod && loda > plodf) {
            curqtl <- qtlspa
            curformula <- formulaspa
            curplod <- ploda
            lod <- loda+lod
            curnqtl <- n.qtl+2
          }
          else if(plodf > curplod) {
            curqtl <- qtlspf
            curformula <- formulaspf
            curplod <- plodf
            lod <- lodf+lod
            curnqtl <- n.qtl+2
          }
        }
        
      }
      
      qtl <- curqtl
      n.qtl <- curnqtl
      attr(qtl, "formula") <- deparseQTLformula(curformula)
      attr(qtl, "pLOD") <- curplod
      formula <- curformula
      lod <- curlod
      
      if(refine.locations) {
        if(verbose) cat(" ---Refining positions\n")
        rqtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl,
                          covar=covar, formula=formula, method=method,
                          verbose=verbose.scan, incl.markers=incl.markers,
                          keeplodprofile=FALSE)
        if(any(rqtl$pos != qtl$pos)) { # updated positions
          if(verbose) cat(" ---  Moved a bit\n")
          qtl <- rqtl
          lod <- fitqtl(cross, pheno.col, qtl, covar=covar, formula=formula,
                        method=method, model=model, dropone=FALSE, get.ests=FALSE,
                        run.checks=FALSE, tol=tol, maxit=maxit)$result.full[1,4] - lod0
          curplod <- calc.plodQuoc(lod, formula=formula, qtl=qtl,
                               penalties=penalties)
          attr(qtl, "pLOD") <- curplod
        }
        
      }
      
      if(verbose)
        cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
            deparseQTLformula(formula), "\n")
      if(verbose > 1) 
        cat("         qtl:", paste(qtl$chr, round(qtl$pos,1), sep="@"), "\n")
      
      
      if(curplod > curbestplod) {
        if(verbose) 
          cat("** new best ** (pLOD increased by ", round(curplod - curbestplod, 4),
              ")\n", sep="")
        
        curbest <- qtl
        curbestplod <- curplod
      }
      
      if(keeptrace) {
        temp <- list(chr=qtl$chr, pos=qtl$pos)
        attr(temp, "formula") <- deparseQTLformula(formula)
        attr(temp, "pLOD") <- curplod
        class(temp) <- c("compactqtl", "list")
        temp <- list(temp)
        names(temp) <- i
        thetrace <- c(thetrace, temp)
      }
      
      if(n.qtl >= max.qtl) break
    }
    
    if(verbose) cat(" -Starting backward deletion\n")
    
    while(n.qtl > 1) {
      i <- i+1
      out <- fitqtl(cross, pheno.col, qtl, covar=covar, formula=formula,
                    method=method, model=model, dropone=TRUE, get.ests=FALSE,
                    run.checks=FALSE, tol=tol, maxit=maxit)$result.drop 
      
      rn <- rownames(out)
      # ignore things with covariates
      wh <- c(grep("^[Qq][0-9]+$", rn),
              grep("^[Qq][0-9]+:[Qq][0-9]+$", rn))
      out <- out[wh,,drop=FALSE]
      
      thelod <- out[,3]
      minlod <- min(thelod, na.rm=TRUE)
      
      wh <- which(!is.na(thelod) & thelod==minlod)
      if(length(wh) > 1) wh <- sample(wh, 1)
      
      lod <- lod - minlod
      todrop <- rownames(out)[wh]
      
      if(verbose) cat(" ---Dropping", todrop, "\n")
      
      if(length(grep(":", todrop)) > 0) { # dropping an interaction
        theterms <- attr(terms(formula), "factors")
        wh <- colnames(theterms)==todrop
        if(!any(wh)) stop("Confusion about what interation to drop!")
        theterms <- colnames(theterms)[!wh]
        formula <- as.formula(paste("y~", paste(theterms, collapse="+")))
      }
      else {
        numtodrop <- as.numeric(substr(todrop, 2, nchar(todrop)))
        
        theterms <- attr(terms(formula), "factors")
        cn <- colnames(theterms)
        g <- c(grep(paste("^[Qq]", numtodrop, "$", sep=""), cn),
               grep(paste("^[Qq]", numtodrop, ":", sep=""), cn),
               grep(paste(":[Qq]", numtodrop, "$", sep=""), cn))
        cn <- cn[-g]
        formula <- as.formula(paste("y~", paste(cn, collapse="+")))
        
        if(n.qtl > numtodrop) {
          for(j in (numtodrop+1):n.qtl)
            formula <- reviseqtlnuminformula(formula, j, j-1)
        }
        
        qtl <- dropfromqtl(qtl, index=numtodrop)
        qtl$name <- qtl$altname <- paste("Q", 1:qtl$n.qtl, sep="")
        n.qtl <- n.qtl - 1
      }
      
      curplod <- calc.plodQuoc(lod, formula=formula, qtl=qtl,
                           penalties=penalties)
      
      if(verbose)
        cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
            deparseQTLformula(formula), "\n")
      if(verbose > 1) 
        cat("         qtl:", paste(qtl$chr, round(qtl$pos,1), sep=":"), "\n")
      
      attr(qtl, "formula") <- deparseQTLformula(formula)
      attr(qtl, "pLOD") <- curplod
      
      if(refine.locations) {
        if(verbose) cat(" ---Refining positions\n")
        if(!is.null(qtl)) {
          rqtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl,
                            covar=covar, formula=formula, method=method,
                            verbose=verbose.scan, incl.markers=incl.markers,
                            keeplodprofile=FALSE)
          if(any(rqtl$pos != qtl$pos)) { # updated positions
            if(verbose) cat(" ---  Moved a bit\n")
            qtl <- rqtl
            lod <- fitqtl(cross, pheno.col, qtl, covar=covar, formula=formula,
                          method=method, model=model, dropone=FALSE, get.ests=FALSE,
                          run.checks=FALSE, tol=tol, maxit=maxit)$result.full[1,4] - lod0
            curplod <- calc.plodQuoc(lod, formula=formula, qtl=qtl,
                                 penalties=penalties)
            attr(qtl, "pLOD") <- curplod
          }
        }
      }
      
      if(curplod > curbestplod) {
        if(verbose) 
          cat("** new best ** (pLOD increased by ", round(curplod - curbestplod, 4),
              ")\n", sep="")
        
        curbestplod <- curplod
        curbest <- qtl
      }
      
      if(keeptrace) {
        temp <- list(chr=qtl$chr, pos=qtl$pos)
        attr(temp, "formula") <- deparseQTLformula(formula)
        attr(temp, "pLOD") <- curplod
        class(temp) <- c("compactqtl", "list")
        temp <- list(temp)
        names(temp) <- i
        thetrace <- c(thetrace, temp)
      }
    }
    
    # re-form the qtl
    if(!is.null(curbest)) {
      chr <- curbest$chr
      pos <- curbest$pos
      o <- order(factor(chr, levels=names(cross$geno)), pos)
      
      qtl <- makeqtl(cross, chr[o], pos[o], what=qtlmethod)
      
      # need to redo numbering in formula
      formula <- as.formula(attr(curbest, "formula"))
      
      if(length(chr) > 1) {
        n.qtl <- length(chr)
        for(i in 1:n.qtl)
          formula <- reviseqtlnuminformula(formula, i, n.qtl+i)
        for(i in 1:n.qtl)
          formula <- reviseqtlnuminformula(formula, n.qtl+o[i], i)
      }
      
      if(keeplodprofile) {
        if(verbose) cat(" ---One last pass through refineqtl\n")
        qtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl,
                         covar=covar, formula=formula, method=method,
                         verbose=verbose.scan, incl.markers=incl.markers,
                         keeplodprofile=TRUE)
      }
      attr(qtl, "formula") <- deparseQTLformula(formula)
      attr(qtl, "pLOD") <- attr(curbest, "pLOD")
      curbest <- qtl
    }
    else {
      curbest <- numeric(0)
      class(curbest) <- "qtl"
      attr(curbest,"pLOD") <- 0
    }
    
    if(keeptrace)
      attr(curbest, "trace") <- thetrace
    
    attr(curbest, "formula") <- deparseQTLformula(attr(curbest, "formula"), TRUE)
    
    curbest
  }

######################################################################
# penalized LOD score, Quoc changed, now central function for pLOD
######################################################################
calc.plodQuoc <-
  function(lod, nterms, type=c("f2","bc"), penalties, formula, qtl) {
    if (missing(nterms) & !missing(formula) & !missing(qtl)) 
      nterms <- countqtltermsQuoc(formula, qtl)
    nterms <- nterms[1:6]
    if(any(penalties==Inf & nterms > 0)) return(-Inf)
    
    as.numeric(lod - sum((nterms*penalties)[nterms > 0]))
  }

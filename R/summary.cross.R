######################################################################
#
# summary.cross.R
#
# copyright (c) 2001-8, Karl W Broman
# last modified Jun, 2008
# first written Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: summary.cross, print.summary.cross, nind, nchr, nmar,
#           totmar, nphe, nmissing, print.cross, chrlen
#
######################################################################

summary.cross <-
function(object,...)
{
  if(!any(class(object) == "cross"))
    stop("Input should have class \"cross\".")
    
  n.ind <- nind(object)
  tot.mar <- totmar(object)
  n.phe <- nphe(object)
  n.chr <- nchr(object)
  n.mar <- nmar(object)
  type <- class(object)[1]

  if(type != "f2" && type != "bc" && type != "4way" &&
     type != "riself" && type != "risib" && type != "cc" && type != "dh") 
    stop("Cross type ", type, " is not suppoted.")

  if(type=="cc") {
    cat("A Collaborative cross\n")
    return(NULL)
  }

  # combine genotype data into one big matrix
  Geno <- pull.geno(object)

  # proportion of missing genotype data
  missing.gen <- mean(is.na(Geno))
  
  # table of genotype values
  if(type=="f2") {
    typings <- table(factor(Geno[!is.na(Geno)], levels=1:5))
    temp <- getgenonames("f2", "A", cross.attr=attributes(object))
    names(typings) <- c(temp, paste("not", temp[c(3,1)]))
  }
  else if(type=="bc" || type=="riself" || type=="risib" || type=="dh") {
    typings <- table(factor(Geno[!is.na(Geno)], levels=1:2))
    names(typings) <- getgenonames(type, "A", cross.attr=attributes(object))
  }
  else if(type=="4way") {
    typings <- table(factor(Geno[!is.na(Geno)], levels=1:14))

    temp <- getgenonames("4way", "A", cross.attr=attributes(object))
    names(typings) <- c(temp,
                        paste(temp[c(1,3)], collapse="/"),
                        paste(temp[c(2,4)], collapse="/"),
                        paste(temp[c(1,2)], collapse="/"),
                        paste(temp[c(3,4)], collapse="/"),
                        paste(temp[c(1,4)], collapse="/"),
                        paste(temp[c(2,3)], collapse="/"),
                        paste("not", temp[1]),
                        paste("not", temp[2]),
                        paste("not", temp[3]),
                        paste("not", temp[4]))
  }
  else typings <- table(factor(Geno[!is.na(Geno)]))

  # turn into fractions
  typings <- typings/sum(typings)

  # amount of missing phenotype data
  missing.phe <- as.numeric(cbind(apply(object$pheno,2,function(a) mean(is.na(a)))))

  # check that, in the case of a "4way" cross, the genetic
  #     maps are matrices with 2 rows, and that for other crosses,
  #     the genetic maps are numeric vectors
  if(type=="4way") {
    if(any(!sapply(object$geno, function(a) (is.matrix(a$map) && nrow(a$map)==2)))) 
      warning("The genetic maps should all be matrices with two rows.")
  }
  else {
    if(any(sapply(object$geno, function(a) is.matrix(a$map))))
      warning("The genetic maps should all be numeric vectors rather than matrices.")
  }

  # check that object$geno[[i]]$data has colnames and that they match
  #     the names in object$geno[[i]]$map
  jitterwarning <- 0
  for(i in 1:n.chr) {
    nam1 <- colnames(object$geno[[i]]$data)
    map <- object$geno[[i]]$map
    if(is.matrix(map)) nam2 <- colnames(map)
    else nam2 <- names(map)
    chr <- names(object$geno)[[i]]
    if(is.null(nam1)) {
      warn <- paste("The data matrix for chr", chr,
                    "lacks column names")
      warning(warn)
    }
    if(is.null(nam2)) {
      warn <- paste("The genetic map for chr", chr,
                    "lacks column names")
      warning(warn)
    }
    if(any(nam1 != nam2)) 
      stop("Marker names in the data matrix and genetic map\n",
                    "for chr ", chr, " do not match.")
      
    if((is.matrix(map) && (any(diff(map[1,])<0) || any(diff(map[2,])<0))) ||
       (!is.matrix(map) && any(diff(map)<0))) 
      stop("Markers out of order on chr ", chr)


    # check that no two markers are on top of each other
    if(is.matrix(map)) { # sex-specific maps
      n <- ncol(map)
      if(n > 1) {
        d1 <- diff(map[1,])
        d2 <- diff(map[2,])
        if(any(d1 < 1e-14 & d2 < 1e-14))
          jitterwarning <- 1
      }
    }
    else {
      n <- length(map)
      if(n > 1) {
        d <- diff(map)
        if(any(d < 1e-14)) 
          jitterwarning <- 1
      }
    }

  }
    
  if(jitterwarning)
    warning("Some markers at the same position; use jittermap().")


  if(!is.data.frame(object$pheno)) 
    warning("Phenotypes should be a data.frame.")

  x <- table(colnames(object$pheno))
  if(any(x > 1)) 
    warning("Some phenotypes have the same name:\n",
            paste(names(x)[x>1], collapse="  "))

  # check genotype data
  if(type=="bc" || type=="riself" || type=="risib" || type=="dh") {
    # Invalid genotypes?
    if(any(!is.na(Geno) & Geno != 1 & Geno != 2)) { 
      u <- unique(as.numeric(Geno))
      u <- sort(u[!is.na(u)])
      warn <- paste("Invalid genotypes.",
                    "\n    Observed genotypes:",
                    paste(u, collapse=" "))
      warning(warn)
    }

    # Missing genotype category on autosomes?
    if(sum(!is.na(Geno) & Geno==2) == 0 ||
       sum(!is.na(Geno) & Geno==1) == 0) {
      warn <- paste("Strange genotype pattern on chr ", chr, ".", sep="")
      warning(warn)
    }
  }
  else if(type=="f2") {
    # invalid genotypes
    if(any(!is.na(Geno) & Geno!=1 & Geno!=2 & Geno!=3 &
           Geno!=4 & Geno!=5)) { 
      u <- unique(as.numeric(Geno))
      u <- sort(u[!is.na(u)])
      warn <- paste("Invalid genotypes on chr", chr, ".", 
                    "\n    Observed genotypes:",
                    paste(u, collapse=" "))
      warning(warn)
    }

    # X chromosome
    for(i in 1:n.chr) {
      if(class(object$geno[[i]]) == "X") {
        dat <- object$geno[[i]]$data
        if(any(!is.na(dat) & dat!=1 & dat!=2)) { 
          u <- unique(as.numeric(dat))
          u <- sort(u[!is.na(u)])
          warn <- paste("Invalid genotypes on X chromosome:",
                        "\n    Observed genotypes:",
                        paste(u, collapse=" "))
          warning(warn)
        }
      }
    }

    # Missing genotype category on autosomes?
    dat <- NULL; flag <- 0
    for(i in 1:n.chr) {
      if(class(object$geno[[i]]) != "X") {
        dat <- cbind(dat,object$geno[[i]]$data)
        flag <- 1
      }
    }
    if(flag && (sum(!is.na(dat) & dat==2) == 0 ||
               sum(!is.na(dat) & dat==1) == 0 ||
               sum(!is.na(dat) & dat==3) == 0))
      warning("Strange genotype pattern.")
  }
  else if(type=="4way") {
    # Invalid genotypes?
    if(any(!is.na(Geno) & (Geno != round(Geno) | Geno < 1 | Geno > 14))) {
      u <- unique(as.numeric(Geno))
      u <- sort(u[!is.na(u)])
      warn <- paste("Invalid genotypes.",
                    "\n    Observed genotypes:",
                    paste(u, collapse=" "))
      warning(warn)
    }
  }


  # Look for duplicate marker names
  mnames <- NULL
  for(i in 1:nchr(object)) 
    mnames <- c(mnames,colnames(object$geno[[i]]$data))

  o <- table(mnames)
  if(any(o > 1))
    warning("Duplicate markers [", paste(names(o)[o>1], collapse=", "), "]")

  # make sure the genotype data are matrices rather than data frames
  if(any(sapply(object$geno, function(a) is.data.frame(a$data))))
    warning("The $data objects should be simple matrices, not data frames.")

  # make sure each chromosome has class "A" or "X"
  chr.class <- sapply(object$geno, class)
  if(!all(chr.class == "A" | chr.class == "X"))
    warning("Each chromosome should have class \"A\" or \"X\".")
  chr.nam <- names(object$geno)
  if(any(chr.class=="A" & (chr.nam=="X" | chr.nam=="x"))) {
    wh <- which(chr.nam=="X" | chr.nam=="x")
    warning("Chromosome \"", chr.nam[wh], "\" has class \"A\" but probably ",
            "should have class \"X\".")
  }

  # if more than one X chromosome, print a warning
  if(sum(chr.class=="X") > 1)
    warning("More than one X chromosome: [",
            paste(chr.nam[chr.class == "X"], collapse=", "), "]")

  if(any(chr.class=="A")) 
    autosomes <- chr.nam[chr.class == "A"]
  else autosomes <- NULL
  if(any(chr.class=="X"))
    Xchr <- chr.nam[chr.class == "X"]
  else Xchr <- NULL

  # check individual IDs
  id <- getid(object)
  if(!is.null(id) && length(id) != length(unique(id)))
    warning("The individual IDs are not unique.")

  cross.summary <- list(type=type, n.ind = n.ind, n.phe=n.phe, 
			n.chr=n.chr, n.mar=n.mar,
			missing.gen=missing.gen,typing.freq=typings,
			missing.phe=missing.phe,
                        autosomes=autosomes, Xchr=Xchr)
  class(cross.summary) <- "summary.cross"
  cross.summary
  
}


print.summary.cross <-
function(x,...)
{
#  cat("\n")
  if(x$type=="f2") cat("    F2 intercross\n\n")
  else if(x$type=="bc") cat("    Backcross\n\n")
  else if(x$type=="4way") cat("    4-way cross\n\n")
  else if(x$type=="riself") cat("    RI strains via selfing\n\n")
  else if(x$type=="risib") cat("    RI strains via sib matings\n\n")
  else if(x$type=="dh") cat("    Doubled haploids\n\n")
  else cat(paste("    cross", x$type, "\n\n",sep=" "))

  cat("    No. individuals:   ", x$n.ind,"\n\n")
  cat("    No. phenotypes:    ", x$n.phe,"\n")
  cat("    Percent phenotyped:", round((1-x$missing.phe)*100,1), "\n\n")
  cat("    No. chromosomes:   ", x$n.chr,"\n")
  if(!is.null(x$autosomes))
    cat("        Autosomes:     ", paste(x$autosomes, collapse=" "), "\n")
  if(!is.null(x$Xchr))
    cat("        X chr:         ", paste(x$Xchr, collapse=" "), "\n")
  cat("\n")
  cat("    Total markers:     ", sum(x$n.mar), "\n")
  cat("    No. markers:       ", x$n.mar, "\n")
  cat("    Percent genotyped: ", round((1-x$missing.gen)*100,1), "\n")
  cat("    Genotypes (%):     ", 
      paste(names(x$typing.freq),round(x$typing.freq*100,1),sep=":", collapse="  "),
      "\n")
#  cat("\n")
}



nind <-
function(object)
{
  if(!any(class(object) == "cross"))
    stop("Input should have class \"cross\".")

  n.ind1 <- nrow(object$pheno)
  n.ind2 <- sapply(object$geno,function(x) nrow(x$data))
  if(any(n.ind2 != n.ind1))
    stop("Different numbers of individuals in genotypes and phenotypes.")
  n.ind1
}

nchr <-
function(object)
{
  if(!any(class(object) == "cross"))
    stop("Input should have class \"cross\".")

  length(object$geno)
}

nmar <- 
function(object)
{
  if(!any(class(object) == "cross"))
    stop("Input should have class \"cross\".")

  if(!is.matrix(object$geno[[1]]$map))
    n.mar1 <- sapply(object$geno, function(x) length(x$map))
  else # sex-specific maps
    n.mar1 <- sapply(object$geno, function(x) ncol(x$map))

  n.mar2 <- sapply(object$geno, function(x) ncol(x$data))
  if(any(n.mar1 != n.mar2))
    stop("Different numbers of markers in genotypes and maps.")
  n.mar1
}

totmar <-
function(object)
{
  if(!any(class(object) == "cross"))
    stop("Input should have class \"cross\".")

  if(!is.matrix(object$geno[[1]]$map))
    totmar1 <- sum(sapply(object$geno, function(x) length(x$map)))
  else # sex-specific maps
    totmar1 <- sum(sapply(object$geno, function(x) ncol(x$map)))
  totmar2 <- sum(sapply(object$geno, function(x) ncol(x$data)))
  if(totmar1 != totmar2)
    stop("Different numbers of markers in genotypes and maps.")
  totmar1
}

nphe <-
function(object)
{
  if(!any(class(object) == "cross"))
    stop("Input should have class \"cross\".")

  ncol(object$pheno)
}

# count number of missing genotypes for each individual or each marker
nmissing <-
function(cross,what=c("ind","mar"))
{
  if(!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")

  what <- match.arg(what)

  if(what=="ind") {
    n.missing <- rep(0,nind(cross))
    for(i in 1:nchr(cross)) 
      n.missing <- n.missing +
        apply(cross$geno[[i]]$data,1,function(a) sum(is.na(a)))

    # individual IDs
    id <- getid(cross)
    if(!is.null(id)) names(n.missing) <- id
  }
  else {
    n.missing <- NULL
    for(i in 1:nchr(cross))
      n.missing <- c(n.missing,
                     apply(cross$geno[[i]]$data,2,function(a) sum(is.na(a))))
  }

  n.missing
}


# "print" method for cross object
#
# to avoid ever printing the entire object, print just a little
#     warning message and then the summary
print.cross <-
function(x, ...)
{
  cat("  This is an object of class \"cross\".\n")
  cat("  It is too complex to print, so we provide just this summary.\n")
  print(summary(x))
  return(summary(x))
}


# get chromosome lengths
chrlen <-
function(object)
{
  if(!any(class(object) == "map") && !any(class(object) == "cross"))
    stop("Input should have class \"map\" or \"cross\".")
  
  if(!any(class(object) == "map"))
    x <- pull.map(object)
  else x <- object

  if(is.matrix(x[[1]])) 
    return(sapply(x, apply, 1, function(a) diff(range(a))))

  sapply(x, function(a) diff(range(a)))
}


# end of summary.cross.R

#####################################################################
#
# mqmcircleplot.R
# mainroutine: mqmplot_circle
#
# copyright (c) 2009, Danny Arends
# last modified December, 2009
# first written December, 2009
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# This file contains the following routines:
# 
#   test_mqmplot_circle     : Test / initialization routine to check the circle functions / needed libraries+datasets
#   mqmplot_circle          : Circleplotting routine for scanone / mqmmulti
#   circlelocations         : Helperfunction to subdivide a circle
#   drawspline              : Draws a spline between 2 locations on the circle
#   getchromosomelength     : Retrieve chromosome length from scanone object
#   getgenomelength         : Retrieve genome length form scanone object
#   locationtocircle        : Translates a position in cM to a Circular position
#   drawcirculargenome      : Drawing routine for the genome uses a scanon object
#   loopthroughmulti        : Loop through an mqmmultiobject highlighting traits 1 by 1
#
######################################################################

test_mqmplot_circle <- function(){
  library(qtl)
  data(hyper)
  data(multitrait)
  data(locations)
  
  hyperfilled <- fill.geno(hyper)
  hypercof <- mqmsetcofactors(hyper,4)
  hyperres <- mqm(hyperfilled,hypercof)

  multifilled <- fill.geno(multitrait)
  multicof <- mqmsetcofactors(multitrait,10)
  multiloc <- addloctocross(multifilled,locations)
  multires <- mqmall(multifilled,cofactors=multicof)

  mqmplot_circle(hyperfilled,hyperres)
  Sys.sleep(3)
  mqmplot_circle(multifilled,multires)
  Sys.sleep(3)
  mqmplot_circle(multiloc,multires)
  Sys.sleep(3)
  mqmplot_circle(multitrait,multires,highlight=16)
  Sys.sleep(3)
  mqmplot_circle(multiloc,multires,highlight=16)
}

test_plot_circle <- function(){
  values <- (runif(10)*3)+1
  plot.circle(1:10)
  plot.circle(values)
  plot.circle(values,dist=T)
  values <- (runif(1000)*3)+1
  plot.circle(values,dist=T,cex=0.9)
}


plot.circle <- function(x, distance=F, vulcano=F, 
                        point.color="blue", text.color="red", line.color=rgb(0.1,0.1,0.1,0.1), 
                        lwd=1, cex=1, spacing=1, verbose=FALSE){
  if(is.null(x)) stop("Please supply a vector in x")
  if(is.vector(x)){
    ploty <- c(-1.1,1.1)
    if(vulcano) ploty <- c(0.0,3.0)
    plot(c(-1.1, 1.1), ploty, type = "n", axes = F, xlab = "", ylab = "")
    title(main = "Circular plot")
    cvalues <- circlelocations(length(x)+spacing*length(x))
    l <- 1
    itemnum <- 1
    for(item in x){
      nl <- l + spacing
      loc <- c(0,0)
      if(vulcano) loc <- c(0,1.75)
      if(!distance){
        loc <- loc+t(cvalues[l+(nl-l),])
      }else{
        if(is.numeric(x)){
          loc <- loc+(item/max(x)) * t(cvalues[l+(nl-l),])
        }
        if(is.character(x)){
          loc <- loc+(nchar(item)/max(nchar(x))) * t(cvalues[l+(nl-l),])
        }
        points(loc,cex=0.5*cex,col=point.color)
      }
      drawspline(c(0,0),loc,col=line.color,lwd=lwd)
      if(!is.na(text.color)){
        if(is.null(names(x))){
          text(1.1*loc,paste("",format(item,dig=2)),col=text.color,cex=0.7*cex)
        }else{
          text(1.1*loc,paste(names(x)[itemnum]),col=text.color,cex=0.7*cex)
        }
      }
      l <- 1 + nl
      itemnum <- itemnum + 1
    }
  }else{
    stop("Please supply a vector in x")
  }
}

mqmplotcircle <- function(...){
  mqmplot_circle(...)
}

mqmplot_circle <- function(cross,result,highlight=0,spacing=25,legend=FALSE,verbose=FALSE){
  if(is.null(cross)){
		stop("No cross object. Please supply a valid cross object.") 
	}
  lod<-FALSE
  if(any(class(result)=="mqmmulti")){
    if(highlight > 0){
      templateresult <- mqmextractmarkers(result[[highlight]])
      lod<-TRUE
    }else{
      templateresult <- mqmextractmarkers(result[[1]])
      lod<-FALSE
    }
  }else{
    templateresult <- mqmextractmarkers(result)
    lod<-TRUE
  }
  if(!("scanone" %in% class(templateresult))){
    stop("Wrong type of result file, please supply a valid scanone object.") 
  }
  colorz <- rainbow(length(result),alpha=0.8)
  if(!legend){
    totallength <- getgenomelength(templateresult)
    nchr <- length(unique(templateresult[,1]))
    cvalues <- circlelocations(totallength+(nchr*spacing))

    drawcirculargenome(templateresult,spacing=spacing,lod=lod)
    
    if(any(class(result)=="mqmmulti")){
      #multiple scan results, so plot em all
      for(x in 1:length(result)){
        model <- mqmgetmodel(result[[x]])
        if(!is.null(model)){
          if(verbose) cat("Model found for trait ",x,"\n")
          if(!is.null(cross$locations)){
            if(verbose) cat("Locations of traits available\n")
            title(sub = "Multiple traits with locations")
            location <- as.numeric(cross$locations[[x]])
            traitl <- locationtocircle(templateresult,location[1],location[2],spacing=spacing)
          }else{
            if(verbose) cat("Locations of traits not available\n")
            title(sub = "Multiple traits no locations")
            traitl <- t(c(0,0+0.25*(0.5-(x/length(result)))))
          }
          if(highlight==0){
            col=colorz[x]
          }else{
            col=rgb(0.1, 0.1, 0.1, 0.1)  
          }
          points(traitl,col=col,pch=24,cex=1)
          for(y in 1:length(model[[4]])){
            qtll <- locationtocircle(templateresult,model[[4]][y],model[[5]][y],spacing=spacing)
            drawspline(traitl,qtll,col=col)
            points(qtll*(1+0.1*((x/length(result)))),col=col,pch=19,cex=1)
          }
        }else{
          if(verbose) cat("Trait ",x," has no model\n")
        }
      }
      legend("bottomleft",c("Trait","QTL"),col=c("black","black"),pch=c(24,19),cex=1)
    }
    if(any(class(result)=="scanone") || highlight > 0){
      #single scan result or highlighting one of the multiple
      if(!any(class(result)=="scanone")){
        #Just highlight the template result
        result <- templateresult
      }
      if(verbose) cat("Single trait plot\n")
      model <- mqmgetmodel(result)
      if(!is.null(model)){
         for(y in 1:length(model[[4]])){
          qtll <- locationtocircle(templateresult,model[[4]][y],model[[5]][y],spacing=spacing)
          points(qtll,col="red",pch=19,cex=1)
          text(qtll*1.15,model[[2]][y],col="red",cex=0.7)
          if(!is.null(cross$locations)){
            location <- as.numeric(cross$locations[[1]])
            traitl <- locationtocircle(templateresult,location[1],location[2],spacing=spacing)
            points(traitl,col="red",pch=24,cex=1)
          }else{
            traitl <- c(0,0)
          }
          drawspline(traitl,qtll,col="red")
        }   
      }
      legend("bottomright","Significant Cofactor",col="red",pch=19,cex=0.75)
      if(highlight==0) title(sub = "Single trait")
    }
  }else{
    plot(c(-1,1), c(-1, 1), type = "n", axes = F, xlab = "", ylab = "")
    title(main = "Legend to circular genome plot")
    legend("center",paste(colnames(cross$pheno)),col=colorz,pch=19,cex=0.75)
  }
}

circlelocations <- function(nt){
  medpoints <- matrix(nrow = nt, ncol = 2)
  phi <- seq(0, 2 * pi, length = (nt + 1))
  complex.circle <- complex(modulus = 1, argument = phi)
  for (j in 1:nt) {
    medpoints[j, ] <- c(Im(complex.circle[j]), Re(complex.circle[j]))
  }
  medpoints
}

drawspline <- function (cn1, cn2, lwd = 1,col="blue",...){
  x <- cbind(cn1[1],0,cn2[1])
  y <- cbind(cn1[2],0,cn2[2])
  xspline(x, y, shape=1, lwd=lwd, border=col,...)
}

getchromosomelength <- function(result, chr){
  l <- ceiling(max(result[which(result[,1]==chr),2]))
  l
}

getgenomelength <- function(result){
  l <- 1
  for(x in unique(result[,1])){
    l <- l + getchromosomelength(result,x)
  }
  l
}

locationtocircle <- function(result, chr, loc, spacing=50, fixoutofbounds=TRUE, verbose=FALSE){
  templateresult <- result
  totallength <- getgenomelength(result)
  nchr <- length(unique(templateresult[,1]))
  cvalues <- circlelocations(totallength+(nchr*spacing))
  l <- 1
  for(x in unique(templateresult[,1])){
    if(x==chr){
      if(loc < getchromosomelength(result,x)){
        return(t(cvalues[(l+loc),]))
      }else{
        if(verbose) cat("Location out of chromosome bounds",loc," ",getchromosomelength(result,x),"\n")
        if(fixoutofbounds) return(t(cvalues[(l+getchromosomelength(result,x)),]))
        stop(paste("Location out of chromosome bounds",loc," ",getchromosomelength(result,x),"\n"))
      }
    }
    l <- l + getchromosomelength(result,x) + spacing
  }
  stop("No such chromosome")
}

drawcirculargenome <- function(result,lodmarkers=FALSE,spacing=50){
  result <- mqmextractmarkers(result)
  plot(c(-1.1, 1.1), c(-1.1, 1.1), type = "n", axes = F, xlab = "", ylab = "")
  title(main = "Circular genome plot")
  totallength <- getgenomelength(result)
  nchr <- length(unique(result[,1]))
  cvalues <- circlelocations(totallength+(nchr*spacing))
  l <- 1
  for(x in unique(result[,1])){
    #Draw chromosomes
    nl <- l+getchromosomelength(result,x) 
    lines(cvalues[l:nl,],cex=0.01)
    l <- nl + spacing
  }
  for(x in 1:nrow(result)){
    #Draw markers
    if(lodmarkers){
      size <- min(c((result[x,3]/2+1),4))
      c <- gray(0.5-(0.4*(result[x,3]/max(result[,3]))))
      points(locationtocircle(result,result[x,1],result[x,2],spacing=spacing),pch=20,col=c,cex=size)
    }else{
      points(locationtocircle(result,result[x,1],result[x,2],spacing=spacing),pch=20)
    }
  }
}

loopthroughmulti <- function(cross,result,save=FALSE,spacing=100){
  n <- 1
  while(n <= length(multires)){
    if(save) jpeg(file=paste("circleplotT",n,".jpg",sep=""),w=800,h=800)
    mqmplot_circle(cross,result,spacing=spacing,highlight=n)
    if(save) dev.off()
    n <- n+1
  }
}


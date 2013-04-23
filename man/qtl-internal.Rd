\name{qtl-internal}
\alias{calc.plod}
\alias{checkcovar}
\alias{clean}
\alias{condense}
\alias{convert}
\alias{convert2bcsft}
\alias{create.map}
\alias{find.pseudomarkerpos}
\alias{flipcross}
\alias{genAllPartitions}
\alias{getsex}
\alias{addmarkerstointervalmap}
\alias{mqm_version}
\alias{scanall}
\alias{getgenonames}
\alias{imf.cf}
\alias{imf.h}
\alias{imf.k}
\alias{imf.m}
\alias{imf.stahl}
\alias{locatemarker}
\alias{max.scantwocondensed}
\alias{mf.cf}
\alias{mf.h}
\alias{mf.k}
\alias{mf.m}
\alias{mf.stahl}
\alias{polyplot}
\alias{print.addcovarint}
\alias{print.addint}
\alias{print.compactqtl}
\alias{print.cross}
\alias{print.map}
\alias{print.qtl}
\alias{print.scanoneboot}
\alias{print.scantwo}
\alias{print.summary.addpair}
\alias{print.summary.compactqtl}
\alias{print.summary.cross}
\alias{print.summary.fitqtl}
\alias{print.summary.map}
\alias{print.summary.qtl}
\alias{print.summary.ripple}
\alias{print.summary.scanone}
\alias{print.summary.scanoneperm}
\alias{print.summary.scantwo}
\alias{print.summary.scantwo.old}
\alias{print.summary.scantwoperm}
\alias{replacemap}
\alias{subrousummaryscantwo}
\alias{summary.scantwocondensed}
\alias{summary.addcovarint}
\alias{summary.addint}
\alias{summary.compactqtl}

\title{Internal qtl functions}

\description{
  Internal qtl functions. These are generally not to be called by the
  user.
}

\usage{
calc.plod(lod, nterms, type=c("f2","bc"), penalties)
checkcovar(cross, pheno.col, addcovar, intcovar, perm.strata,
           ind.noqtl=NULL, weights=NULL, verbose=TRUE)
clean(object, \dots)
condense(object)
convert(object, \dots)
convert2bcsft(cross, BC.gen = 0, F.gen = 0, estimate.map = TRUE,
              error.prob=0.0001, map.function=c("haldane","kosambi","c-f","morgan"),
              \dots)
create.map(map, step, off.end, stepwidth=c("fixed", "variable", "max"))
find.pseudomarkerpos(cross, marker, where=c("draws", "prob"))
flipcross(cross)
genAllPartitions(n.taxa, taxa)
getsex(cross)
addmarkerstointervalmap(cross,intervalresult,verbose=FALSE)
mqm_version()
scanall(cross, scanfunction=scanone, multicore=TRUE, n.clusters=1, batchsize=10, FF=0,cofactors=NULL, \dots, plot=FALSE, verbose=FALSE)
getgenonames(type=c("f2","bc","riself","risib","4way","dh","special","bcsft"),
             chrtype=c("A","X"), expandX=c("simple","standard","full"),
             sexpgm, cross.attr)
imf.cf(r)
imf.h(r)
imf.k(r)
imf.m(r)
imf.stahl(r, m=0, p=0, tol=1e-12, maxit=1000)
locatemarker(map, pos, chr, flag)
\method{max}{scantwocondensed}(object, lodcolumn=1,
    what=c("best", "full", "add", "int"),
    df=FALSE, na.rm=TRUE, \dots)
mf.cf(d)
mf.h(d)
mf.k(d)
mf.m(d)
mf.stahl(d, m=0, p=0)
polyplot(x, type='b', legend=TRUE,legendloc=0, labels=NULL, cex =
         par("cex"), pch = 19, gpch = 21, bg = par("bg"), color = par("fg"),
         col=NULL, ylim=range(x[is.finite(x)]), xlim = NULL, main = NULL, xlab =
         NULL, ylab = NULL, add=FALSE, \dots) 
\method{print}{addcovarint}(x, \dots)
\method{print}{addint}(x, \dots)
\method{print}{compactqtl}(x, \dots)
\method{print}{cross}(x, \dots)
\method{print}{map}(x, \dots)
\method{print}{qtl}(x, \dots)
\method{print}{scanoneboot}(x, \dots)
\method{print}{scantwo}(x, \dots)
\method{print}{summary.addpair}(x, \dots)
\method{print}{summary.compactqtl}(x, \dots)
\method{print}{summary.cross}(x, \dots)
\method{print}{summary.fitqtl}(x, \dots)
\method{print}{summary.map}(x, \dots)
\method{print}{summary.qtl}(x, \dots)
\method{print}{summary.ripple}(x, \dots)
\method{print}{summary.scanone}(x, \dots)
\method{print}{summary.scanoneperm}(x, \dots)
\method{print}{summary.scantwo}(x, \dots)
\method{print}{summary.scantwo.old}(x, \dots)
\method{print}{summary.scantwoperm}(x, \dots)
replacemap(object, map)
subrousummaryscantwo(object, for.perm=FALSE)
\method{summary}{scantwocondensed}(object, thresholds, 
         what=c("best", "full", "add", "int"),
         perms, alphas, lodcolumn=1, pvalues=FALSE,
         df=FALSE, allpairs=TRUE, \dots)
\method{summary}{addcovarint}(object, \dots)
\method{summary}{addint}(object, \dots)
\method{summary}{compactqtl}(object, \dots)
}


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{internal}

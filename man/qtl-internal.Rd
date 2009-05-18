\name{qtl-internal}
\alias{adjust.rf.ri}
\alias{asnumericwithdec}
\alias{calc.genoprob.special}
\alias{calc.pairprob}
\alias{calc.plod}
\alias{charround}
\alias{checkcovar}
\alias{checkformula}
\alias{checkStepwiseqtlStart}
\alias{clean}
\alias{condense}
\alias{convert}
\alias{convertMWril}
\alias{countqtlterms}
\alias{create.map}
\alias{deparseQTLformula}
\alias{discan}
\alias{dropfromqtlformula}
\alias{dropXcol}
\alias{effectplot.getmark}
\alias{effectplot.calmeanse}
\alias{expandf2covar}
\alias{find.pseudomarkerpos}
\alias{fitqtlengine}
\alias{fitstahl.estp}
\alias{fitstahl.este}
\alias{fitstahl.estpe}
\alias{fitstahl.estp.sub}
\alias{fitstahl.este.sub}
\alias{fitstahl.estpe.sub}
\alias{fixX4write}
\alias{fixXgeno.bc}
\alias{fixXgeno.f2}
\alias{genotab.em}
\alias{getid}
\alias{getsex}
\alias{getgenonames}
\alias{grab.arg.names}
\alias{imf.cf}
\alias{imf.h}
\alias{imf.k}
\alias{imf.m}
\alias{imf.stahl}
\alias{interpmap}
\alias{interpmap4scantwo}
\alias{LikePheVector}
\alias{locatemarker}
\alias{locate.xo}
\alias{makeSSmap}
\alias{markerforwsel}
\alias{markerforwself2}
\alias{markerloglik}
\alias{matchchr}
\alias{max.scantwocondensed}
\alias{mf.cf}
\alias{mf.h}
\alias{mf.k}
\alias{mf.m}
\alias{mf.stahl}
\alias{mybinaryrep}
\alias{orderMarkers.sub}
\alias{parseformula}
\alias{plot.effectscan}
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
\alias{printQTLformulanicely}
\alias{qtlformulasymmetric}
\alias{read.cro.qtlcart}
\alias{read.cross.csv}
\alias{read.cross.csvs}
\alias{read.cross.gary}
\alias{read.cross.karl}
\alias{read.cross.mm}
\alias{read.cross.qtlcart}
\alias{read.cross.qtx}
\alias{read.map.qtlcart}
\alias{read.maps.mm}
\alias{reorgRIargmax}
\alias{reorgRIdraws}
\alias{reorgRIgenoprob}
\alias{reorgRIpairprob}
\alias{replacemap}
\alias{revisecovar}
\alias{reviseqtlnuminformula}
\alias{revisescantwodf}
\alias{reviseXdata}
\alias{ripple.perm1}
\alias{ripple.perm2}
\alias{ripple.perm.sub}
\alias{roundqtlpos}
\alias{scanone.perm}
\alias{scanone.perm.engine}
\alias{scanoneXnull}
\alias{scantwo.perm}
\alias{scantwo.perm.engine}
\alias{sim.bcg}
\alias{sim.cross.4way}
\alias{sim.cross.bc}
\alias{sim.cross.f2}
\alias{sim.ril}
\alias{subrousummaryscantwo}
\alias{summary.scantwocondensed}
\alias{summary.addcovarint}
\alias{summary.addint}
\alias{summary.compactqtl}
\alias{testchr}
\alias{vbscan}
\alias{write.cross.csv}
\alias{write.cross.gary}
\alias{write.cross.mm}
\alias{write.cross.qtlcart}

\title{Internal qtl functions}

\description{
  Internal qtl functions. These are generally not to be called by the
  user.
}

\usage{
adjust.rf.ri(r, type=c("self","sib"), chrtype=c("A","X"), expand=TRUE)
asnumericwithdec(x, dec=".")
calc.genoprob.special(cross, error.prob=0.0001,
                      map.function=c("haldane","kosambi","c-f","morgan"))
calc.pairprob(cross, step=0, off.end=0, error.prob=0.0001, 
              map.function=c("haldane","kosambi","c-f","morgan"), map,
              assumeCondIndep=FALSE)
calc.plod(lod, nterms, type=c("f2","bc"), penalties)
charround(x, digits=1)
checkcovar(cross, pheno.col, addcovar, intcovar, perm.strata, verbose=TRUE)
checkformula(formula, qtl.name, covar.name)
checkStepwiseqtlStart(qtl, formula, covar=NULL)
clean(object, \dots)
condense(object)
convert(object)
convertMWril(cross, parents, error.prob=0)
countqtlterms(formula, ignore.covar=TRUE)
create.map(map, step, off.end, stepwidth=c("fixed", "variable"))
deparseQTLformula(formula, reorderterms=FALSE)
discan(cross, pheno, method=c("em","mr"),
       addcovar=NULL, intcovar=NULL, maxit=4000, tol=1e-4,
       verbose=FALSE, give.warnings=TRUE)
dropfromqtlformula(formula, qtlnum)
dropXcol(type=c("f2","bc"), sexpgm, cross.attr)
effectplot.getmark(cross, mname)
effectplot.calmeanse(pheno, mark1, mark2, geno1, geno2, ndraws, var.flag=c("pooled", "group"))
expandf2covar(thecovar)
find.pseudomarkerpos(cross, marker, where=c("draws", "prob"))
fitqtlengine(pheno, qtl, covar=NULL, formula, method=c("imp", "hk"),
             dropone=TRUE, get.ests=FALSE, run.checks=TRUE, cross.attr)
fitstahl.estp(cross, error.prob=0.0001, m=0, tol=1e-4, maxit=4000)
fitstahl.estp.sub(p, cross, error.prob=0.0001, m=0, thetol=1e-4, maxit=4000)
fitstahl.este(cross, m=0, p=0, tol=1e-4, maxit=4000)
fitstahl.este.sub(error.prob, cross, m=0, p=0, thetol=1e-4, maxit=4000)
fitstahl.estpe(cross, m=0, tol=1e-4, maxit=4000)
fitstahl.estpe.sub(x, cross, m=0, thetol=1e-4, maxit=4000)
fixX4write(geno,sex,pgm,crosstype)
fixXgeno.bc(cross)
fixXgeno.f2(cross, alleles)
markerforwsel(x, y, maxsize=7)
markerforwself2(x, y, maxsize=7)
genotab.em(dat, tol=1e-6, maxit=10000, verbose=FALSE)
getid(cross)
getsex(cross)
getgenonames(type=c("f2","bc","riself","risib","4way","dh","special"),
             chrtype=c("A","X"), expandX=c("simple","standard","full"),
             sexpgm, cross.attr)
grab.arg.names(\dots)
imf.cf(r)
imf.h(r)
imf.k(r)
imf.m(r)
imf.stahl(r, m=0, p=0, tol=1e-12, maxit=1000)
interpmap(oldmap, newmap)
interpmap4scantwo(output, newmap)
LikePheVector(pheno, n.ind, n.phe)
locatemarker(map, pos, chr, flag)
locate.xo(cross)
makeSSmap(cross)
markerloglik(cross, marker, error.prob=0.0001)
matchchr(selection, thechr)
\method{max}{scantwocondensed}(object, lodcolumn=1,
    what=c("best", "full", "add", "int"),
    df=FALSE, na.rm=TRUE, \dots)
mf.cf(d)
mf.h(d)
mf.k(d)
mf.m(d)
mf.stahl(d, m=0, p=0)
mybinaryrep(n)
orderMarkers.sub(cross, chr, window=7, use.ripple=TRUE)
parseformula(formula, qtl.dimname, covar.dimname)
\method{plot}{effectscan}(x, gap=25, ylim, mtick=c("line","triangle"),
     add.legend=TRUE, alternate.chrid=FALSE, \dots)
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
printQTLformulanicely(formula, header, width, width2, sep=" ")
qtlformulasymmetric(formula, qtlnum1, qtlnum2)
read.cro.qtlcart(file)
read.cross.csv(dir, file, na.strings=c("-","NA"),
               genotypes=c("A","H","B","D","C"),
               estimate.map=TRUE, rotate=FALSE, \dots)
read.cross.csvs(dir, genfile, phefile, na.strings=c("-","NA"),
                genotypes=c("A","H","B","D","C"),
                estimate.map=TRUE, rotate=FALSE, \dots)
read.cross.gary(dir, genfile, mnamesfile, chridfile, phefile,
                pnamesfile, mapfile,estimate.map,na.strings)
read.cross.karl(dir, genfile, mapfile, phefile)
read.cross.mm(dir, rawfile, mapfile, estimate.map=TRUE)
read.cross.qtlcart(dir, crofile, mapfile)
read.cross.qtx(dir, file, estimate.map=TRUE)
read.map.qtlcart(file)
read.maps.mm(mapsfile)
reorgRIargmax(cross)
reorgRIdraws(cross)
reorgRIgenoprob(cross)
reorgRIpairprob(cross, pairprob)
replacemap(object, map)
revisecovar(sexpgm, covar)
reviseqtlnuminformula(formula, oldnum, newnum)
revisescantwodf(df)
reviseXdata(type=c("f2","bc"), expandX=c("simple","standard","full"),
            sexpgm, geno, prob, draws, pairprob, cross.attr)
ripple.perm1(n)
ripple.perm2(n)
ripple.perm.sub(x, mat)
roundqtlpos(x, digits = 1)
scanone.perm(cross, pheno.col=1, model=c("normal","binary","2part","np"),
             method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use=c("all.obs", "complete.obs"), upper=FALSE,
             ties.random=FALSE, start=NULL, maxit=4000,
             tol=1e-4, n.perm=1000, perm.Xsp=FALSE, perm.strata=NULL,
             verbose=TRUE, batchsize=250)
scanone.perm.engine(n.perm, cross, pheno.col, model,
                    method, addcovar, intcovar, weights, use,
                    upper, ties.random, start, maxit, tol,
                    verbose, perm.strata, batchsize=250)
scanoneXnull(type, sexpgm)
scantwo.perm(cross, pheno.col=1, model=c("normal","binary"),
             method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use=c("all.obs", "complete.obs"),
             incl.markers=FALSE, clean.output=FALSE, 
             clean.nmar=1, clean.distance=0,
             maxit=4000, tol=1e-4, verbose=FALSE,
             n.perm=1000, perm.strata, assumeCondIndep=FALSE,
             batchsize=250)
scantwo.perm.engine(n.perm, cross, pheno.col, model,
                    method, addcovar, intcovar, weights, use,
                    incl.markers, clean.output,  
                    clean.nmar=1, clean.distance=0,
                    maxit, tol, verbose,
                    perm.strata, assumeCondIndep=FALSE,
                    batchsize=250)
sim.bcg(n.ind, map, m, p,
        map.function=c("haldane","kosambi","c-f","morgan"))
sim.cross.4way(map, model, n.ind, error.prob, missing.prob,
               partial.missing.prob, keep.errorind, m, p, map.function)
sim.cross.bc(map, model, n.ind, error.prob, missing.prob,
             keep.errorind, m, p, map.function)
sim.cross.f2(map, model, n.ind, error.prob, missing.prob,
             partial.missing.prob, keep.errorind, m, p, map.function)
sim.ril(map, n.ril=1, type=c("sibmating", "selfing"), n.str=c("2","4","8"),
        m=0, p=0, error.prob=0, missing.prob=0, random.cross=TRUE)
subrousummaryscantwo(object, for.perm=FALSE)
\method{summary}{scantwocondensed}(object, thresholds, 
         what=c("best", "full", "add", "int"),
         perms, alphas, lodcolumn=1, pvalues=FALSE,
         df=FALSE, allpairs=TRUE, \dots)
\method{summary}{addcovarint}(object, \dots)
\method{summary}{addint}(object, \dots)
\method{summary}{compactqtl}(object, \dots)
testchr(selection, thechr)
vbscan(cross, pheno.col=1, upper=FALSE, method="em", maxit=4000,
       tol=1e-4)
write.cross.csv(cross, filestem="data", digits=5, rotate=FALSE,
                split=FALSE)
write.cross.gary(cross, digits)
write.cross.mm(cross, filestem="data", digits=5)
write.cross.qtlcart(cross, filestem="data")
}


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{internal}

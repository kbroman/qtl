\name{read.cross}
\alias{read.cross}

\title{Read data for a QTL experiment}

\description{
  Data for a QTL experiment is read from a set of files and converted
  into an object of class \code{cross}.  The comma-delimited format
  (\code{csv}) is recommended.  All formats require chromosome
  assignments for the genetic markers, and assume that markers are in
  their correct order.
}

\usage{
read.cross(format=c("csv", "csvr", "csvs", "csvsr", "mm", "qtx",
                    "qtlcart", "gary", "karl"),
           dir="", file, genfile, mapfile, phefile, chridfile,
           mnamesfile, pnamesfile, na.strings=c("-","NA"),
           genotypes=c("A","H","B","D","C"), alleles=c("A","B"),
           estimate.map=TRUE, convertXdata=TRUE, error.prob=0.0001,
           map.function=c("haldane", "kosambi", "c-f", "morgan"),
           BC.gen=0, F.gen=0, \dots)
}

\arguments{
  \item{format}{Specifies the format of the data.}
  \item{dir}{Directory in which the data files will be found.  In
    Windows, use forward slashes (\code{"/"}) or double backslashes
    (\code{"\\\\"}) to specify directory trees.}
  \item{file}{The main input file for formats \code{csv}, \code{csvr}
    and \code{mm}.} 
  \item{genfile}{File with genotype data (formats \code{csvs},
    \code{csvsr}, \code{karl} and \code{gary} only).}
  \item{mapfile}{File with marker position information (all 
    except the \code{csv} formats).}
  \item{phefile}{File with phenotype data (formats \code{csvs},
    \code{csvsr}, \code{karl} and \code{gary} only).}
  \item{chridfile}{File with chromosome ID for each marker (\code{gary}
    format only).} 
  \item{mnamesfile}{File with marker names (\code{gary} format only).}
  \item{pnamesfile}{File with phenotype names (\code{gary} format
    only).}
  \item{na.strings}{A vector of strings which are to be interpreted as
    missing values (\code{csv} and \code{gary} formats only).  For the
    \code{csv} formats, these are interpreted globally
    for the entire
    file, so missing value codes in phenotypes must not be valid
    genotypes, and vice versa. For the \code{gary} format, these are
    used only for the phenotype data.} 
  \item{genotypes}{A vector of character strings specifying the genotype
    codes (\code{csv} formats only).  Generally this is a vector of
    length 5, with the elements corresponding to AA, AB, BB, not BB 
    (i.e., AA or AB), and not AA (i.e., AB or BB).  \bold{Note}: Pay
    careful attention to the third and fourth of these; the order of
    these can be confusing.

    If you are trying to read 4-way cross data, your file must have
    genotypes coded as described below, and you need to set 
    \code{genotypes=NULL} so that no re-coding gets done.}
  \item{alleles}{A vector of two one-letter character strings (or four,
    for the four-way cross), to be used as labels for the two alleles.}
  \item{estimate.map}{For all formats but \code{qtlcart} and
    \code{karl}: if TRUE and marker positions are
    not included in the input files, the genetic map is estimated using
    the function  \code{\link{est.map}}.} 
  \item{convertXdata}{If TRUE, any X chromosome genotype data is
    converted to the internal standard, using columns \code{sex} and
    \code{pgm} in the phenotype data if they available or by inference
    if they are not.  If FALSE, the X chromsome data is read as is.}  
  \item{error.prob}{In the case that the marker map must be estimated:
   Assumed genotyping error rate used in the calculation of the
   penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{In the case that the marker map must be estimated:
   Indicates whether to use the Haldane, Kosambi,
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions. (Ignored if m > 0.)}
  \item{BC.gen}{Used only for cross type \code{"bcsft"}.}
  \item{F.gen}{Used only for cross type \code{"bcsft"}.}
  \item{\dots}{Additional arguments, passed to the function
    \code{\link[utils]{read.table}} in the case of 
    \code{csv} and \code{csvr} formats.  In particular, one may use the
    argument
    \code{sep} to specify the field separator (the default is a comma),
    \code{dec} to specify the character used for the decimal point
    (the default is a period), and \code{comment.char} to specify a
    character to indicate comment lines.}
}

\value{
  An object of class \code{cross}, which is a list with two components:

  \item{geno}{This is a list with elements corresponding to
    chromosomes.  \code{names(geno)} contains the names of the
    chromsomes.  Each chromosome is itself a list, and is given class
    \code{A} or \code{X} according to whether it is autosomal
    or the X chromosome. 
    
    There are two components for each chromosome: \code{data}, a matrix
    whose rows are individuals and whose columns are markers, and
    \code{map}, either a vector of marker positions (in cM) or a matrix
    of dim (\code{2 x n.mar}) where the rows correspond to marker
    positions in female and male genetic distance, respectively.  

    The genotype data gets converted into numeric codes, as follows.  

    The genotype data for a backcross is coded as NA = missing,
    1 = AA, 2 = AB.  

    For an F2 intercross, the coding is NA = missing, 1 = AA, 2 = AB, 3
    = BB, 4 = not BB (i.e. AA or AB; D in Mapmaker/qtl), 5 = not AA (i.e. AB
    or BB; C in Mapmaker/qtl). 

    For a 4-way cross, the mother and father are assumed to have
    genotypes AB and CD, respectively.  The genotype data for the
    progeny is assumed to be phase-known, with the following coding
    scheme: NA = missing, 1 = AC, 2 = BC, 3 = AD, 4 = BD, 5 = A = AC or AD,
    6 = B = BC or BD, 7 = C = AC or BC, 8 = D = AD or BD, 9 = AC or BD,
    10 = AD or BC, 11 = not AC, 12 = not BC, 13 = not AD, 14 = not BD.
  }

  \item{pheno}{data.frame of size (\code{n.ind x n.phe}) containing the
    phenotypes.  If a phenotype with the name \code{id} or \code{ID} is
  included, these identifiers will be used in \code{\link{top.errorlod}},
  \code{\link{plotErrorlod}}, and \code{\link{plotGeno}} as
  identifiers for the individual.} 

  While the data format is complicated, there are a number of functions,
  such as \code{\link{subset.cross}}, to
  assist in pulling out portions of the data.
}

\details{
  The available formats are comma-delimited (\code{csv}), rotated
  comma-delimited (\code{csvr}), comma-delimited with separate files for
  genotype and phenotype data (\code{csvs}), rotated comma-delimited
  with separate files for genotype and phenotype data (\code{csvsr}),
  Mapmaker (\code{mm}), Map Manager QTX (\code{qtx}), Gary Churchill's
  format (\code{gary}) and Karl Broman's format (\code{karl}).  The
  required files and their specification for each format appears below.
  The comma-delimited formats are recommended. Note that most of these
  formats work only for backcross and intercross data.

  The \code{sampledata} directory in the package distribution contains
  sample data files in multiple formats.    Also see
  \url{http://www.rqtl.org/sampledata}.

  The \code{\dots} argument enables additional arguments to be passed to
  the function \code{\link[utils]{read.table}} in the case of \code{csv}
  and \code{csvr} formats.  In particular, one may use the argument
  \code{sep} to specify the field separator (the default is a comma),
  \code{dec} to specify the character used for the decimal point (the
  default is a period), and \code{comment.char} to specify a character
  to indicate comment lines.
}


\section{X chromosome}{
  \bold{The genotypes for the X chromosome require special care!}

  The X chromosome should be given chromosome identifier \code{X} or
  \code{x}.  If it is labeled by a number or by \code{Xchr}, it will be
  interpreted as an autosome.

  The phenotype data should contain a column named \code{"sex"} which
  indicates the sex of each individual, either coded as \code{0}=female and
  \code{1}=male, or as a factor with levels \code{female}/\code{male} or
  \code{f}/\code{m}.  Case will be 
  ignored both in the name and in the factor levels.  If no such
  phenotype column is included, it will be assumed that all individuals
  are of the same sex.

  In the case of an intercross, the phenotype data may also contain a
  column named \code{"pgm"} (for "paternal grandmother") indicating the
  direction of the cross.  It should be coded as 0/1 with 0 indicating
  the cross (AxB)x(AxB) or (BxA)x(AxB) and 1 indicating the cross
  (AxB)x(BxA) or (BxA)x(BxA).  If no such phenotype column is included,
  it will be assumed that all individuals come from the same direction
  of cross. 

  The internal storage of X chromosome data is quite different from that
  of autosomal data.  Males are coded 1=AA and 2=BB; females with pgm==0
  are coded 1=AA and 2=AB; and females with pgm==1 are coded 1=BB and
  2=AB.  If the argument \code{convertXdata} is TRUE, conversion to this
  format is made automatically; if FALSE, no conversion is done, 
  \code{\link{summary.cross}} will likely return a warning, and
  most analyses will not work properly.

  Use of \code{convertXdata=FALSE} (in which case the X chromosome
  genotypes will not be converted to our internal standard) can be
  useful for diagnosing problems in the data, but will require some
  serious mucking about in the internal data structure.
}

\section{CSV format}{
  The input file is a comma-delimited text file.  A different field
  separator may be specified via the argument \code{sep}, which will be passed
  to the function \code{\link[utils]{read.table}}).  For example, in
  Europe, it is common to use a comma in place of the decimal point in
  numbers and so a semi-colon in place of a comma as the field
  separator; such data may be read by using \code{sep=";"} and
  \code{dec=","}.  

  The first line should contain the phenotype names followed by the
  marker names.  \bold{At least one phenotype must be included}; for
  example, include a numerical index for each individual.

  The second line should contain blanks in the phenotype columns,
  followed by chromosome identifiers for each marker in all other
  columns. If a chromosome has the identifier \code{X} or \code{x}, it
  is assumed to be the X chromosome; otherwise, it is assumed to be an
  autosome.

  An optional third line should contain blanks in the phenotype
  columns, followed by marker positions, in cM.

  Marker order is taken from the cM positions, if provided; otherwise,
  it is taken from the column order.

  Subsequent lines should give the data, with one line for each
  individual, and with phenotypes followed by genotypes.  If possible,
  phenotypes are made numeric; otherwise they are converted to factors.

  The genotype codes must be the same across all markers.  For example,
  you can't have one marker coded AA/AB/BB and another coded A/H/B.
  This includes genotypes for the X chromosome, for which hemizygous
  individuals should be coded as if they were homoyzogous.

  The cross is determined to be a backcross if only the first two elements
  of the \code{genotypes} string are found; otherwise, it is assumed to
  be an intercross.
}

\section{CSVr format}{
  This is just like the \code{csv} format, but rotated (or really
  transposed), so that rows are columns and columns are rows.
}


\section{CSVs format}{
  This is like the \code{csv} format, but with separate files for the
  genotype and phenotype data.

  The first column in the genotype data must specify individuals'
  identifiers, and there must be a column in the phenotype data with
  precisely the same information (and with the same name).  These IDs
  will be included in the data as a phenotype.  If the name \code{id} or
  \code{ID} is used, these identifiers will be used in
  \code{\link{top.errorlod}}, \code{\link{plotErrorlod}}, and
  \code{\link{plotGeno}} as identifiers for the individual.  

  The first row in each file contains the column names.  For the
  phenotype file, these are the names of the phenotypes.  For the
  genotype file, the first cell will be the name of the identifier
  column (\code{id} or \code{ID}) and the subsequent fields will be the
  marker names.

  In the genotype data file, the second row gives the chromosome IDs.
  The cell in the second row, first column, must be blank.  A third
  row giving cM positions of markers may be included, in which case the
  cell in the third row, first column, must be blank.

  There need be no blank rows in the phenotype data file.
}


\section{CSVsr format}{
  This is just like the \code{csvs} format, but with each file rotated
  (or really transposed), so that rows are columns and columns are rows.
}


\section{Mapmaker format}{
  This format requires two files.  The so-called rawfile, specified by
  the argument \code{file}, contains the genotype and phenotype
  data. Rows beginning with the symbol \code{#} are ignored.  The first
  line should be either \code{data type f2 intercross} or
  \code{data type f2 backcross}.  The second line should begin with
  three numbers indicating the numbers of individuals, markers and
  phenotypes in the file.  This line may include the word \code{symbols}
  followed by symbol assignments (see the documentation for mapmaker,
  and cross your fingers).  The rest of the lines give genotype data
  followed by phenotype data, with marker and phenotype names always
  beginning with the \code{*} symbol.
    
  A second file contains the genetic map information, specified with
  the argument \code{mapfile}.  The map file may be in
  one of two formats.  The function will determine which format of map
  file is presented.  

  The simplest format for the map file is not standard for the Mapmaker
  software, but is easy to create.  The file contains two or three
  columns separated by white space and with no header row.  The first
  column gives the chromosome assignments.  The second column gives the
  marker names, with markers listed in the order along the chromosomes.
  An optional third column lists the map positions of the markers.

  Another possible format for the map file is the \code{.maps}
  format, which is produced by Mapmaker.  The code for reading this
  format was written by Brian Yandell.

  Marker order is taken from the map file, either by the order they are
  presented or by the cM positions, if specified.  
}

\section{Map Manager QTX format}{
  This format requires a single file (that produced by the Map Manager
  QTX program).
}  

\section{QTL Cartographer format}{
  This format requires two files: the \code{.cro} and \code{.map} files
  for QTL Cartographer (produced by the QTL Cartographer
  sub-program, Rmap and Rcross).

  Note that the QTL Cartographer cross types are converted as follows:
  RF1 to riself, RF2 to risib, RF0 (doubled haploids) to bc, B1 or B2 to
  bc, RF2 or SF2 to f2.
}  

\section{Gary format}{
  This format requires the six files.  All files have default names, and
  so the file names need not be specified if the default names are used. 

  \code{genfile} (default = \code{"geno.dat"}) contains the genotype
  data.  The file contains one line per individual, with genotypes for
  the set of markers separated by white space.  Missing values are
  coded as 9, and genotypes are coded as 0/1/2 for AA/AB/BB.

  \code{mapfile} (default = \code{"markerpos.txt"}) contains two
  columns with no header row: the marker names in the first column and
  their cM position in the second column.  If marker positions are not
  available, use \code{mapfile=NULL}, and a dummy map will be inserted.
  
  \code{phefile} (default = \code{"pheno.dat"}) contains the phenotype
  data, with one row for each mouse and one column for each phenotype.
  There should be no header row, and missing values are coded as
  \code{"-"}. 

  \code{chridfile} (default = \code{"chrid.dat"}) contains the
  chromosome identifier for each marker.  

  \code{mnamesfile} (default = \code{"mnames.txt"}) contains the marker
  names.

  \code{pnamesfile} (default = \code{"pnames.txt"}) contains the names
  of the phenotypes.  If phenotype names file is not available, use
  \code{pnamesfile=NULL}; arbitrary phenotype names will then be
  assigned. 
}


\section{Karl format}{
  This format requires three files; all files have default names, and so
  need not be specified if the default name is used. 

  \code{genfile} (default = \code{"gen.txt"}) contains the genotype
  data.  The file contains one line per individual, with genotypes
  separated by white space.  Missing values are coded 0; genotypes are
  coded as 1/2/3/4/5 for AA/AB/BB/not BB/not AA.
    
  \code{mapfile} (default = \code{"map.txt"}) contains the map
  information, in the following complicated format: \cr \cr
    \code{n.chr} \cr
    \code{n.mar(1) rf(1,1) rf(1,2) \ldots rf(1,n.mar(1)-1)}\cr
    \code{mar.name(1,1)}\cr 
    \code{mar.name(1,2)}\cr 
    \code{\ldots}\cr 
    \code{mar.name(1,n.mar(1))}\cr 
    \code{n.mar(2)}\cr
    \code{\ldots}\cr 
    \code{etc.} \cr

  \code{phefile} (default = \code{"phe.txt"}) contains a matrix of
  phenotypes, with one individual per line.  The first line in the
  file should give the phenotype names.
}

\examples{
\dontrun{
# CSV format
dat1 <- read.cross("csv", dir="Mydata", file="mydata.csv")

# CSVS format
dat2 <- read.cross("csvs", dir="Mydata", genfile="mydata_gen.csv",
                   phefile="mydata_phe.csv")

# you can read files directly from the internet
datweb <- read.cross("csv", "http://www.rqtl.org/sampledata",
                     "listeria.csv")

# Mapmaker format
dat3 <- read.cross("mm", dir="Mydata", file="mydata.raw",
                   mapfile="mydata.map")

# Map Manager QTX format
dat4 <- read.cross("qtx", dir="Mydata", file="mydata.qtx")

# QTL Cartographer format
dat5 <- read.cross("qtlcart", dir="Mydata", file="qtlcart.cro",
                   mapfile="qtlcart.map")

# Gary format
dat6 <- read.cross("gary", dir="Mydata", genfile="geno.dat",
                   mapfile="markerpos.txt", phefile="pheno.dat",
                   chridfile="chrid.dat", mnamesfile="mnames.txt",
                   pnamesfile="pnames.txt")

# Karl format
dat7 <- read.cross("karl", dir="Mydata", genfile="gen.txt",
                   phefile="phe.txt", mapfile="map.txt")}
} 

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}; Brian S. Yandell}

\references{
  Broman, K. W. and Sen,
  \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. (2009) \emph{A
  guide to QTL mapping with R/qtl.}  Springer.  \url{http://www.rqtl.org/book}
}


\seealso{ \code{\link{subset.cross}}, \code{\link{summary.cross}},
  \code{\link{plot.cross}}, \code{\link{c.cross}}, \code{\link{clean.cross}},
  \code{\link{write.cross}}, \code{\link{sim.cross}}, \code{\link[utils]{read.table}}.
  The \code{sampledata} directory in the package distribution contains
  sample data files in multiple formats.  Also see
  \url{http://www.rqtl.org/sampledata}.  
}

\keyword{IO}

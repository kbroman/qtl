#! /bin/sh

echo "* Starting SWEAVE to generate chapter12 (MQM)"
which R
R --version

R CMD BATCH SweaveIt.R
latex chapter12.tex
# run twice for references
latex chapter12.tex 
dvipdfm chapter12.dvi

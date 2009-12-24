#! /bin/sh

echo "* Starting SWEAVE to generate MQM-tour (MQM)"
which R
R --version

R CMD BATCH SweaveIt.R
latex MQM-tour.tex
# run twice for references
latex MQM-tour.tex 
dvipdfm MQM-tour.dvi

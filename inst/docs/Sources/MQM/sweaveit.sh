#! /bin/sh

echo "* Starting SWEAVE to generate MQM-tour (MQM)"
rm -vf MQM-tour.tex
which R
R --version

R CMD BATCH SweaveIt.R
if [ ! -e MQM-tour.tex ]; then
  cat SweaveIt.Rout 
  exit 2
fi
latex MQM-tour.tex
# run twice for references
latex MQM-tour.tex 
dvipdfm MQM-tour.dvi

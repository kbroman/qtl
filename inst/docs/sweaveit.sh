#! /bin/sh

echo "* Starting SWEAVE to generate chapter12 (MQM)"

R CMD BATCH SweaveIt.R
latex chapter12.tex
dvipdfm chapter12.dvi

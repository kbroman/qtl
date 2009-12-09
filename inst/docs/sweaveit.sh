#! /bin/sh

R CMD BATCH SweaveIt.R
latex chapter12.tex
dvipdfm chapter12.dvi

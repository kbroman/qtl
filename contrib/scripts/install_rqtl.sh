#! /bin/sh
#
#  Usage: install_rqtl.sh [path]
#

R CMD build $1
R CMD INSTALL qtl_*.gz

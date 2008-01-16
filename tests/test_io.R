######################################################################
#
# TestIO/input.R
#
# copyright (c) 2002, Karl W Broman
# last modified Feb, 2002
# first written Feb, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# This file contains code for testing the cross IO in R/qtl.
#
# Needed input files:
#
#    gen.txt, map.txt, phe.txt    [Karl's format]
#    listeria.raw, listeria.map   [mapmaker format]
#    listeria.raw, listeria2.map  [mapmaker format; no marker pos]
#    listeria.csv                 [csv format] 
#    listeria2.csv                [csv format; no marker pos] 
#
######################################################################

library(qtl)

##############################
# Reading
##############################
# Read CSV format
csv <- read.cross("csv", "", "listeria.csv")
csv2 <- read.cross("csv", "", "listeria2.csv", estimate=FALSE)

# Read mapmaker format
mm <- read.cross("mm", "", "listeria.raw", "listeria.map")
mm2 <- read.cross("mm", "", "listeria.raw", "listeria2.map", estimate=FALSE)

##############################
# Writing
##############################
# Write in CSV format
write.cross(csv, "csv", filestem="junk1")
csv3 <- read.cross("csv", "", "junk1.csv", genotypes=c("AA","AB","BB","not BB","not AA"))
comparecrosses(csv, csv3)

# Write in mapmaker format
write.cross(csv, "mm", filestem="junk2")

######################################################################
#
# TestIO/input.R
#
# copyright (c) 2002, Karl W Broman
# last modified Feb, 2002
# first written Feb, 2002
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
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

# Cleanup
unlink("junk1.csv")
unlink("junk2.raw")
unlink("junk2.prep")

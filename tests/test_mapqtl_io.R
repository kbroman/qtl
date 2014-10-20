# test input/output in mapqtl format

library(qtl)
data(fake.4way)

write.cross(fake.4way, "mapqtl", "fake_4way_mapqtl")

x <- read.cross("mapqtl", "", genfile="fake_4way_mapqtl.loc",
                phefile="fake_4way_mapqtl.qua",
                mapfile="fake_4way_mapqtl_female.map")

x <- replace.map(x, pull.map(fake.4way))

comparecrosses(x, fake.4way)

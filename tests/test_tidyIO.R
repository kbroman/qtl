library(qtl)
data(hyper)

# write to tidy format
write.cross(hyper, "tidy", "hyper_tidy")

# read back in
x <- read.cross("tidy", "", genfile="hyper_tidy_gen.csv",
                mapfile="hyper_tidy_map.csv", phefile="hyper_tidy_phe.csv",
                genotypes=c("BB", "BA", "AA"))

# compare results
comparecrosses(x, hyper)

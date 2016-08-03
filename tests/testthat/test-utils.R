context("Utility functions")

test_that(desc = 'is.cross()',
          code = {
            x <- 27599

            expect_false(object = is.cross(x = x))

            class(x) <- 'cross'
            expect_false(object = is.cross(x = x))

            names(x) <- 'pheno'
            expect_false(object = is.cross(x = x))

            x <- list(pheno = 1:5, geno = 1:5)
            class(x) <- 'cross'
            expect_false(object = is.cross(x = x))

            y <- qtl::sim.cross(map = qtl::sim.map())
            y[['pheno']] <- y[['pheno']][-17,]
            expect_false(object = is.cross(x = y))

            y <- qtl::sim.cross(map = qtl::sim.map())
            y[['geno']][[1]][['map']] <- y[['geno']][[1]][['map']][-1]
            expect_false(object = is.cross(x = y))

            z <- qtl::sim.cross(map = qtl::sim.map())
            expect_true(object = is.cross(x = z))
          })


test_that(desc = 'is.f2.cross()',
          code = {
            x <- qtl::sim.cross(map = qtl::sim.map(), type = 'bc')
            expect_false(object = is.f2.cross(x = x))

            y <- qtl::sim.cross(map = qtl::sim.map(sex.sp = TRUE), type = '4way')
            expect_false(object = is.f2.cross(x = y))

            z <- qtl::sim.cross(map = qtl::sim.map(), type = 'f2')
            expect_true(object = is.f2.cross(x = z))
          })

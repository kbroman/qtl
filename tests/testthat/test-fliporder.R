context("flip.order")

test_that("flip.order, when applied twice, should get us back to the same thing", {

    data(hyper)
    hyper <- calc.genoprob(hyper, step=1)
    hyper <- sim.geno(hyper, step=10, n.draws=2)
    hyper <- argmax.geno(hyper, step=1)
    hyper <- calc.errorlod(hyper)

    hyperfl <- flip.order(hyper, chr=c(1, 4, 6, 15))
    summary(hyperfl)

    hyperfl2 <- flip.order(hyperfl, chr=c(1, 4, 6, 15))
    summary(hyperfl2)

    # having flipped twice, should be back to where we were
    # (except starting locations for each chromosome map)
    expect_null(comparecrosses(shiftmap(hyper), shiftmap(hyperfl2)))

})

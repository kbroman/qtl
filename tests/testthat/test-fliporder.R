context("flip.order")

test_that("flip.order, when applied twice, should get us back to the same thing", {

    data(hyper)

    # reduce size
    set.seed(53307443)
    hyper <- hyper[,sample(nind(hyper), 8)]

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

test_that("flip.order for 4-way cross", {

    data(fake.4way)

    # reduce size
    set.seed(36461124)
    fake.4way <- fake.4way[,sample(nind(fake.4way), 8)]

    fake.4way <- calc.genoprob(fake.4way, step=1)

    fake.4way.fl <- flip.order(fake.4way, chr=c(1, 4, 6, 15))
    summary(fake.4way.fl)

    fake.4way.fl2 <- flip.order(fake.4way.fl, chr=c(1, 4, 6, 15))
    summary(fake.4way.fl2)

    # having flipped twice, should be back to where we were
    # (except starting locations for each chromosome map)
    expect_null(comparecrosses(shiftmap(fake.4way), shiftmap(fake.4way.fl2)))

})

context("scantwo perms")

test_that("scantwo and scantwopermhk give same results", {

    data(hyper)
    hyper <- calc.genoprob(hyper[c(18,19,"X"),])

    set.seed(92999298)
    out1 <- scantwo(hyper, method="hk", n.perm=3, verbose=FALSE)

    set.seed(92999298)
    out2 <- scantwopermhk(hyper, n.perm=3, verbose=FALSE)

    expect_equivalent(out1, out2)

    # X-chr-specific permutations
    set.seed(92999298)
    out1 <- scantwo(hyper, method="hk", n.perm=3, perm.Xsp=TRUE, verbose=FALSE)

    set.seed(92999298)
    out2 <- scantwopermhk(hyper, n.perm=3, perm.Xsp=TRUE, verbose=FALSE)

    expect_equivalent(out1, out2)

})

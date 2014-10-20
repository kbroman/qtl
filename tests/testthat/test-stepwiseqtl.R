context("stepwiseqtl")

test_that("stepwiseqtl works with X-chr-specific perms", {

    data(fake.f2)
    fake.f2 <- calc.genoprob(fake.f2)
    out2 <- scantwo(fake.f2, method="hk", verbose=FALSE)

    set.seed(17370120)
    operm1 <- scantwopermhk(fake.f2, n.perm=10, verbose=FALSE)
    set.seed(17370120)
    operm2 <- scantwopermhk(fake.f2, n.perm=10, perm.Xsp=TRUE, verbose=FALSE)

    pen1 <- calc.penalties(operm1)
    pen2 <- calc.penalties(operm2)

    out.sq1 <- stepwiseqtl(fake.f2, max.qtl=4, penalties=pen1, method="hk", verbose=FALSE)
    out.sq2 <- stepwiseqtl(fake.f2, max.qtl=4, penalties=pen2, method="hk", verbose=FALSE)
    out.sq3 <- stepwiseqtl(fake.f2, chr=1:19, max.qtl=4, penalties=pen2, method="hk", verbose=FALSE)

})

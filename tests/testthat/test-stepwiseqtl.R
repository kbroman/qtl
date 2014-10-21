context("stepwiseqtl")

test_that("stepwiseqtl works with X-chr-specific perms", {

    data(fake.f2)
    fake.f2 <- calc.genoprob(fake.f2)

    set.seed(17370120)
    operm1 <- scantwopermhk(fake.f2, n.perm=10, verbose=FALSE)
    set.seed(17370120)
    operm2 <- scantwopermhk(fake.f2, n.perm=10, perm.Xsp=TRUE, verbose=FALSE)

    pen1 <- calc.penalties(operm1)
    pen2 <- calc.penalties(operm2)

    out.sq1 <- stepwiseqtl(fake.f2, max.qtl=4, penalties=pen1, method="hk", verbose=FALSE)
    expect_equal(out.sq1$chr, c("1", "8", "13", "X"))
    expect_equal(out.sq1$pos, c(37.11, 61.20, 24.03, 14.20))

    out.sq2 <- stepwiseqtl(fake.f2, max.qtl=4, penalties=pen2, method="hk", verbose=FALSE)
    expect_equal(out.sq2$chr, c("1", "13", "X"))
    expect_equal(out.sq2$pos, c(37.11, 24.03, 14.20))

    out.sq3 <- stepwiseqtl(fake.f2, chr=1:19, max.qtl=4, penalties=pen2, method="hk", verbose=FALSE)
    expect_equal(out.sq3$chr, c("1", "13"))
    expect_equal(out.sq3$pos, c(37.11, 24.03))

    pen2b <- calc.penalties(operm2, alpha=0.2)
    out.sq2b <- stepwiseqtl(fake.f2, max.qtl=6, penalties=pen2b, method="hk", verbose=FALSE)
    expect_equal(out.sq2b$chr, c("1", "8", "13", "X"))
    expect_equal(out.sq2b$pos, c(37.11, 61.20, 24.03, 14.20))

})

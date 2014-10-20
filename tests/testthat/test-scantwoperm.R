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

test_that("summary.scantwo works with X-chr-specific perms", {

    data(hyper)
    set.seed(23615071)
    hyper <- calc.genoprob(fill.geno(hyper[c(18,19,"X"),])) # selected chr; imputed genotypes
    out2 <- scantwo(hyper, method="hk", verbose=FALSE)

    set.seed(17370120)
    operm1 <- scantwopermhk(hyper, n.perm=100, verbose=FALSE)
    set.seed(17370120)
    operm2 <- scantwopermhk(hyper, n.perm=100, perm.Xsp=TRUE, verbose=FALSE)

    # no significant pairs
    sum1 <- summary(out2, perms=operm1, alpha=0.05)
    sum2 <- summary(out2, perms=operm2, alpha=0.05)
    expect_equal(nrow(sum1), 0)
    expect_equal(nrow(sum2), 0)

    # p-values match expectation; not X-chr-specific
    sum1 <- summary(out2, perms=operm1, pvalues=TRUE)
    lodcol <- grep("^lod", names(sum1))
    expect_equal(lodcol, c(5, 7, 9, 13, 15))
    for(i in 1:5)
        expect_equal(sum1[,lodcol[i]+1], sapply(sum1[,lodcol[i]], function(a) mean(operm1[[i]] >= a)))

    # p-values match expectation; X-chr-specific
    sum2 <- summary(out2, perms=operm2, pvalues=TRUE)
    pairtype <- paste0(ifelse(sum2$chr1=="X", "X", "A"),
                       ifelse(sum2$chr2=="X", "X", "A"))
    pairtype <- match(pairtype, c("AA", "AX", "XX"))
    L <- attr(operm2, "L")
    pow <- sum(L)/L
    lodcol <- grep("^lod", names(sum1))
    expect_equal(lodcol, c(5, 7, 9, 13, 15))
    for(j in 1:nrow(sum2)) {
        for(i in 1:5) {
            lod <- sum2[j,lodcol[i]]
            p <- sum2[j,lodcol[i]+1]
            nominal_p <- mean(operm2[[pairtype[j]]][[i]] >= lod)
            adj_p <- 1 - (1-nominal_p)^pow[pairtype[j]]
            expect_equivalent(p, adj_p)
        }
    }

})

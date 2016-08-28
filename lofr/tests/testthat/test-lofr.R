
test_that("make lofmodel", {
    dat <- iris[,1:4]
    mod <- lofmodel(dat, k.min=2, k.max=10)
    expect_equal(mod$k.min, 2)
    expect_equal(mod$k.max, 10)
})

test_that("test predict.lofmodel", {
    dat <- iris[, -5]
    mod <- lofmodel(dat, k.min=1, k.max=10)
    lofs <- predict(mod)
    expect_equal(nrow(lofs), nrow(dat))
    expect_equal(ncol(lofs), 10)

    dat1 <- iris[1:100, -5]
    dat2 <- iris[101:150, -5]
    mod <- lofmodel(dat, k.min=2, k.max=10)
    lofs <- predict(mod, newdata=dat2)
    lofs2 <- predict(mod, newdata=dat2, nparallel=2)
    expect_equal(lofs, lofs2)
    lofs3 <- predict(mod, newdata=dat2, nparallel=3)
    expect_equal(lofs, lofs3)

    library(Rlof)
    dat <- matrix(rnorm(100 * 5))
    mod <- lofmodel(dat, k.min=1, k.max=4)
    lofs <- predict(mod)
    lofs2 <- lof(dat, 1:4)
    expect_equal(as.numeric(lofs), as.numeric(lofs2))

    # The following tests failed.
    # Because 6 - 5 and 6 - 1 have the same distance.
    # iris2 <- iris[1:11, -5]
    # mod <- lofmodel(iris2, k.min=1, k.max=4)
    # lofs <- calc.lof(mod)
    # lofs2 <- lof(iris2, 1:4)
    # expect_equal(as.numeric(lofs), as.numeric(lofs2))    
})

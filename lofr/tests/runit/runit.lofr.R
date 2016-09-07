# -*- coding:utf-8 -*-

test.lofmodel <- function(){
    dat <- iris[,1:4]
    mod <- lofmodel(dat, k.min=2, k.max=10)
    checkEquals(mod$k.min, 2)
    checkEquals(mod$k.max, 10)

    mod2 <- lofmodel(dat, k.min=2, k.max=10, nparallel=2)
    checkEquals(mod$nn.res$nn.index, mod2$nn.res$nn.index)
    checkEquals(mod$nn.res$nn.dist, mod2$nn.res$nn.dist)
    checkEquals(mod$knn.args, mod2$knn.args)
}

test.predict.lofmodel <- function(){
    dat <- iris[, -5]
    mod <- lofmodel(dat, k.min=1, k.max=10)
    lofs <- predict(mod)
    checkEquals(nrow(lofs), nrow(dat))
    checkEquals(ncol(lofs), 10)

    dat1 <- iris[1:100, -5]
    dat2 <- iris[101:150, -5]
    mod <- lofmodel(dat, k.min=2, k.max=10)
    lofs <- predict(mod, newdata=dat2)
    lofs2 <- predict(mod, newdata=dat2, nparallel=2)
    checkEquals(lofs, lofs2)
    lofs3 <- predict(mod, newdata=dat2, nparallel=3)
    checkEquals(lofs, lofs3)

    library(Rlof)
    # Rlof と結果の同一性をみる
    dat <- matrix(rnorm(100 * 5))
    mod <- lofmodel(dat, k.min=1, k.max=4)
    lofs <- predict(mod)
    lofs2 <- lof(dat, 1:4)
    checkEquals(as.numeric(lofs), as.numeric(lofs2))

    # 以下の結果は一致しない
    # k近傍で等距離にある点が複数ある場合の処理に対応していないので
    # 6 - 5 と 6 - 1 が等距離にある。
    # iris2 <- iris[1:11, -5]
    # mod <- lofmodel(iris2, k.min=1, k.max=4)
    # lofs <- calc.lof(mod)
    # lofs2 <- lof(iris2, 1:4)
    # checkEquals(as.numeric(lofs), as.numeric(lofs2))
}

# lofr: Local Outlier Factor for R

## Usage

``` r
library(lofr)
index <- sample(1:150, 100)
dat1  <- iris[index, -5]
dat2  <- iris[setdiff(1:nrow(iris), index), -5]
mod   <- lofmodel(dat1, k.min=1, k.max=10)
lof1  <- predict(mod) # local outier factor for dat1
lof2  <- predict(mod, dat2) # local outlier factor for dat2
```

## Installation

``` r
# install.packages("devtools")
devtools::install_github("5hun/lofr/lofr")
```

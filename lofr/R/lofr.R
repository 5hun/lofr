# -*- coding:utf-8 -*-

#' lofr: A R package for Local Outlier Factor
#'
#' @docType package
#' @name lofr
#' @importFrom FNN get.knn get.knnx
#' @importFrom parallel makeCluster clusterExport parLapply stopCluster
#' @importFrom methods is
NULL

# library(FNN)
# library(parallel)

.LOFMODEL.CLSNAME <- "lofmodel"

.split.data <- function(data, nsplit){
  split.factor <- ceiling(seq_len(nrow(data)) / (nrow(data) / nsplit))
  ret <- list()
  for(j in unique(split.factor)){
    ret <- c(ret, list(data[split.factor == j, , drop=FALSE]))
  }
  return(ret)
}

.get.knn.parallel <- function(data, k, nparallel, knn.args=list()){
  splitted <- .split.data(data, nparallel)
  cl <- makeCluster(nparallel)
  results <- tryCatch({
    clusterExport(cl, "get.knnx")
    results <- parLapply(cl, splitted, 
        function(query, data, k, knn.args){
          knn.args$data  = data
          knn.args$query = query
          knn.args$k     = k
          return(do.call(get.knnx, knn.args))
        }, data=data, k=k+1, knn.args=knn.args)
    results
  }, finally = {
    stopCluster(cl)
  })
  ret <- list(
    nn.index = do.call(rbind, 
        lapply(results, function(x) { x$nn.index[, 2:(k+1), drop=FALSE] })),
    nn.dist  = do.call(rbind,
        lapply(results, function(x) { x$nn.dist[,  2:(k+1), drop=FALSE] }))
  )
  return(ret)
}

.calc.lrd <- function(nn.res, k.min, k.max, nn.res.mod=NULL){
  stopifnot(k.min <= k.max)
  stopifnot(k.min >= 1)
  stopifnot(k.max < nrow(nn.res$nn.index))
  stopifnot(all(dim(nn.res$nn.index) == dim(nn.res$nn.dist)))
  stopifnot(ncol(nn.res$nn.index) == k.max)
  if(is.null(nn.res.mod)){
    nn.res.mod <- nn.res
  }
  stopifnot(all(dim(nn.res.mod$nn.index) == dim(nn.res.mod$nn.dist)))
  stopifnot(ncol(nn.res.mod$nn.index) == k.max)
  stopifnot(all(max(nn.res$nn.index) >= nrow(nn.res.mod$nn.index)))
  
  lrd <- matrix(0.0, nrow=nrow(nn.res$nn.index), ncol=(k.max - k.min + 1))
  for(cur.idx in seq_len(ncol(lrd))){
    cur.k <- cur.idx + k.min - 1
    cur.k.dist <- nn.res.mod$nn.dist[, cur.k]
    k.dist.mat <- matrix(cur.k.dist[c(nn.res$nn.index[, seq_len(cur.k)])],
      nrow=nrow(nn.res$nn.index))
    stopifnot(all(dim(k.dist.mat) == c(nrow(nn.res$nn.index), cur.k)))

    reach.dist <- pmax(k.dist.mat, nn.res$nn.dist[, seq_len(cur.k)])
    lrd[, cur.idx] <- 1 / rowMeans(reach.dist)
  }
  lrd[!is.finite(lrd)] <- 1
  return(lrd)
}

#' Make a model object for \code{\link{predict.lofmodel}}.
#'
#' @param data a matrix.
#' @param k.min mininum value of k.
#' @param k.max maximal value of k.
#' @param knn.args optinal arguments for \code{\link{get.knn}}.
#' @param nparallel optional, the number of cores to be used for parallel computing.
#'
#' @return a list with class "lodmodel".
#'
#' @seealso \code{\link{predict.lofmodel}}, \code{\link[FNN]{get.knn}}.
#'
#' @examples
#' index <- sample(1:150, 100)
#' dat1 <- iris[index, -5]
#' dat2 <- iris[setdiff(1:nrow(iris), index), -5]
#' mod <- lofmodel(dat1, k.min=1, k.max=10)
#' lof1 <- predict(mod) # local outier factor for dat1
#' lof2 <- predict(mod, dat2) # local outlier factor for dat2
#'
#' @export
lofmodel <- function(data, k.min=1, k.max=10, knn.args=list(),
     nparallel=1){
  stopifnot(k.min <= k.max)
  stopifnot(k.min >= 1)
  stopifnot(k.max < nrow(data))

  if(nparallel == 1){
    knn.args$data = data
    knn.args$k = k.max
    nn.res <- do.call(get.knn, knn.args)
  }else{
    nn.res <- .get.knn.parallel(data, k.max, nparallel, knn.args)
  }
    
  stopifnot(all(dim(nn.res$nn.index) == c(nrow(data), k.max)))
  stopifnot(all(dim(nn.res$nn.dist) == c(nrow(data), k.max)))

  lrd <- .calc.lrd(nn.res, k.min, k.max)
  
  model <- list(data=data, k.min=k.min, k.max=k.max, nn.res=nn.res, lrd=lrd, knn.args=knn.args)
  class(model) <- c(.LOFMODEL.CLSNAME, class(model))
  return(model)
}

.calc.lof.from.lrd <- function(lrd.self, nn.res, model){
  k.min <- model$k.min
  k.max <- model$k.max
  lrd.data <- model$lrd

  stopifnot(k.max < nrow(nn.res$nn.indexs))
  stopifnot(ncol(lrd.self) == k.max - k.min + 1)
  stopifnot(nrow(lrd.self) == nrow(nn.res$nn.index))

  lofs <- matrix(0.0, nrow=nrow(lrd.self), ncol=(k.max - k.min + 1))
  for(cur.idx in seq_len(ncol(lofs))){
    cur.k <- cur.idx + k.min - 1
    cur.lrd.self <- lrd.self[, cur.idx]
    cur.lrd.data <- lrd.data[, cur.idx]
    lofs[,cur.idx] <- rowMeans(matrix(cur.lrd.data[nn.res$nn.index[, seq_len(cur.k)]],
        nrow=nrow(lrd.self))) / cur.lrd.self
  }
  return(lofs)
}

.calc.lof.self <- function(model){
  stopifnot(is(model, .LOFMODEL.CLSNAME))
  lof <- .calc.lof.from.lrd(model$lrd, model$nn.res, model)
  return(lof)
}

.calc.lof.other <- function(model, query){
  stopifnot(is(model, .LOFMODEL.CLSNAME))
  stopifnot(ncol(query) == ncol(model$data))

  knn.args <- model$knn.args
  knn.args$query <- query
  nn.res <- do.call(get.knnx, knn.args)
  stopifnot(all(dim(nn.res$nn.index) == c(nrow(query), model$k.max)))
  stopifnot(all(dim(nn.res$nn.dist) == c(nrow(query), model$k.max)))

  lrd <- .calc.lrd(nn.res, model$k.min, model$k.max, model$nn.res)
  lof <- .calc.lof.from.lrd(lrd, nn.res, model)
  return(lof)  
}

#' Calculate local outlier factor.
#'
#' @param object a value of "\code{lofmodel}".
#' @param newdata optional, a matrix to calculate LOF. If omitted, the data given at "\code{lofmodel}" is used.
#' @param nparallel optional, the number of cores to be used for parallel computing.
#' @param ... Additional optional arguments. At present no optional arguments are used.
#'
#' @return a matrix with the local outlier factor of each observation as rows and each k value as columns
#'
#' @details this function does not consider "ties"
#' 
#' @seealso \code{\link{lofmodel}}
#'
#' @examples
#' index <- sample(1:150, 100)
#' dat1 <- iris[index, -5]
#' dat2 <- iris[setdiff(1:nrow(iris), index), -5]
#' mod <- lofmodel(dat1, k.min=1, k.max=10)
#' lof1 <- predict(mod) # local outier factor for dat1
#' lof2 <- predict(mod, dat2) # local outlier factor for dat2
#'
#' @export
predict.lofmodel <- function(object, newdata=NULL, nparallel=1, ...){
  stopifnot(is(object, .LOFMODEL.CLSNAME))
  if(is.null(newdata)){
    stopifnot(nparallel == 1)
    return(.calc.lof.self(object))
  }
  if(nparallel == 1){
    return(.calc.lof.other(object, newdata))
  }else{
    split.factor <- ceiling(seq_len(nrow(newdata)) / (nrow(newdata) / nparallel))
    splitted <- split(newdata, split.factor, drop=FALSE)
    cl <- makeCluster(nparallel)
    clusterExport(cl, 
        c(".LOFMODEL.CLSNAME", "get.knnx", ".calc.lrd", ".calc.lof.from.lrd"))
    results <- parLapply(cl, splitted, .calc.lof.other, model=object)
    stopCluster(cl)
    return(do.call(rbind, results))
  }
}

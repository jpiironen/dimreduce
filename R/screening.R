



# good name: uniscore
#' @export
uniscore <- function(x, y, type='corr', exclude=NULL, normalize=TRUE) {
  
  if (is.vector(x))
    x <- matrix(x, ncol=1)
  
  if (is.factor(y)) {
    # handle factors in a special way
    classes <- levels(y)
    if (length(classes)==2)
      # transform to numeric vector and continue normally
      y <- as.numeric(y==classes[[1]])
    else {
      # more than two class, so score is the maximum from 'class k or other class' scores
      score <- matrix(0, nrow=ncol(x), ncol=1)
      for (c in classes) {
        score <- pmax(score, 
                      uniscore(x, y==c, type=type, exclude=exclude, normalize=normalize))
      }
      return(score)
    }
  }
  
  if (is.vector(y))
    y <- matrix(y, ncol=1)
  
  ok <- rep(TRUE, ncol(x))
  ok[exclude] <- FALSE
  
  if (normalize) {
    x <- scale(x)
    y <- scale(y)
  }
  
  score <- matrix(0, nrow=ncol(x), ncol=ncol(y))
  
  if (type == 'corr') {
    score[ok,] <- abs( t(x[,ok,drop=F]) %*% y )  / (nrow(y)-1)
  } else if (type == 'corr-rank') {
    #apply(x,2,rank)
    stop('not implemented yet.')
  } else if (type == 'runs') {
    score[ok,] <- runs.score(x[,ok,drop=F], y)
  } else
    stop('Unknown score function.')
  
  # if nans produced, treat them as zero score
  score[is.na(score)] <- 0
  
  return(score)
}


# good name: uniscore.test
#' @export
uniscore.test <- function(x,y, type='corr', exclude=NULL, perms=1000, test.max=FALSE) {
  
  # permutations for y
  yperm <- matrix(0, nrow=length(y), ncol=perms)
  for (i in 1:perms)
    yperm[,i] <- sample(y)
  
  # the actual score
  s_actual <- uniscore(x, y, exclude=exclude, type=type)
  
  # compute the scores with permuted y
  if (is.factor(y)) {
    # factors need to be handled separately to make the computation efficient
    classes <- levels(y)
    s_perm <- 0
    for (c in classes) {
      s_perm <- pmax(uniscore(x, yperm==c, exclude=exclude, type=type),
                     s_perm)
    }
  } else {
    s_perm <- uniscore(x,yperm,exclude=exclude, type=type)
  }
  
  if (test.max) {
    pval <- mean( apply(s_perm,2,'max') >= max(s_actual)  )
    return(pval)
  }
  
  eps <- 1e-12 # for numerical stability when there are ties
  pval <- rowMeans(s_perm >= as.vector(s_actual)-eps)
  return(pval)
}





runs.score <- function(x, y) {
  
  if (is.factor(y))
    y <- as.numeric(y)
  if (is.vector(y))
    y <- matrix(y)
  if (is.vector(x))
    x <- matrix(x)
  n <- nrow(y)
  
  if (!is.integer(y)) {
    # bin the values of y into two; to those above median and to those below
    large <- apply(y, 2, function(yj) yj>median(yj))
    y[large] <- 1
    y[!large] <- 0
  }
  
  count <- matrix(0, nrow=ncol(x), ncol=ncol(y)) 
  for (j in 1:ncol(x)) {
    ord <- order(x[,j], runif(n))
    ytemp <- y[ord,,drop=F]
    count[j,] <- colSums(ytemp[1:(n-1),,drop=F] == ytemp[2:n,,drop=F])
  }
  count <- count/(n-1)
  
  if (ncol(count) == 1 || nrow(count) == 1)
    count <- as.vector(count)
  return(count)
}


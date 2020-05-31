
#' Iterative supervised principal components
#'
#' Computes dimension reduction based on the iterative supervised principal components
#' algorithm.
#' 
#' @param x The original feature matrix, columns denoting the features and rows the instances.
#' @param y A vector with the observed target values we try to predict using \code{x}.
#' Can be factor for classification problems.
#' @param nctot Total number of latent features to extract.
#' @param ncsup Maximum number of latent features to extract that use supervision.
#' @param exclude Columns (variables) in x to ignore when extrating the new features.
#' @param nthresh Number of evaluations when finding the optimal screening threshold at each supervised
#' iteration. Increasing this number can make the supervised iterations more accurate but also increases
#' the computation time.
#' @param thresh Instead of specifying \code{nthresh}, one can specify the candidate screening thresholds
#' explicitly. These are numbers between 0 and 1 and are relative to the highest univariate score.
#' By default seq(0,0.999,len=nthresh).
#' @param window Maximum number of features to consider when computing each supervised component.
#' Lowering this number makes the computation faster, but can make the algorithm less accurate
#' if there are more potentially relevant features than this number.
#' @param verbose Whether to print some messages along the way. 
#' @param min_score Terminate the computation at the latest when the maximum univariate score drops
#'  below this.
#' @param normalize Whether to scale the extracted features so that they all have standard deviation
#'  of one.
#' @param center Whether to center the original features before the computation.
#' @param scale Whether to scale the original features to have unit variance before the computation.
#' @param permtest Whether to use permutation test to decide the number of supervised components.
#' @param permtest_type Either 'max-marginal' or 'marginal'.
#' @param alpha Significance level used in the permutation test to decide whether to continue
#' supervised iteration.
#' @param perms Number of permutations to estimate the p-values for univariate scores.
#' @param method Method to compute the principal components. Either 'svd' (default), 
#' 'power' or 'robust'. 'power' can sometimes
#' be slightly faster than 'svd' but in some cases can have very slow convergence. 'robust' 
#' uses robust PCA instead of the standard PCA (requires package pcaPP).
#' @param ... Currently ignored.
#'
#' @return ispca-object that is similar in spirit to the object returned by \code{\link[stats]{prcomp}}.
#' The object will have the following elements:
#' \describe{
#'  \item{\code{w}}{The projection (or rotation) matrix W, that transforms the original data 
#'  \eqn{X} into the new features \eqn{Z = X W} .}
#'  \item{\code{z}}{The extracted latent features corresponding to the training inputs \eqn{X}.}
#'  \item{\code{v}}{Matrix \eqn{V} that is used to compute \eqn{W}. The columns of \eqn{V} indicate
#'  which variables become active at each iteration (see the paper below for more information).}
#'  \item{\code{sdev}}{Standard deviations of the new features.}
#'  \item{\code{ncsup}}{How many supervised components were extracted (the rest are computed 
#'  in an unsupervised manner).}
#'  \item{\code{centers}}{Mean values for the original variables.}
#'  \item{\code{scales}}{Scales of the original variables.}
#'  \item{\code{exclude}}{Excluded variables.}
#' }
#' 
#' @section References:
#' 
#' Piironen, J. and Vehtari, A. (2018). Iterative supervised principal components.
#' In \emph{Proceedings of the 21st International Conference on Artificial
#' Intelligence and Statistics (AISTATS) PMLR 84: 106-114}.
#'
#' @examples
#' \donttest{
#' ### 
#' dr <- ispca(x,y, nctot=2)
#' z <- predict(dr, x) # the latent features
#' }
#'

#' @export
ispca <- function(x,y, nctot=NULL, ncsup=NULL, exclude=NULL, nthresh=NULL, thresh=NULL,
                  window=500, verbose=TRUE, min_score=1e-4, normalize=FALSE,
                  center=TRUE, scale=TRUE, permtest=TRUE, permtest_type='max-marginal', 
                  alpha=0.1, perms=500, method='svd', ...) {
  
  
  if (is.factor(y) && permtest_type!='marginal' && permtest_type!='max-marginal' )
    stop('Requested permutation test for factor type y not implemented yet.')
  
  
  n <- NROW(x)
  d <- NCOL(x)
  
  if (is.null(nctot)) {
    nctot <- min(n-1,d)
    compute_unsup <- FALSE
  } else {
    compute_unsup <- TRUE
  }
  if (nctot > min(n-1,d))
    stop('nctot cannot exceed min(n-1,d)')
  if (is.null(ncsup) || ncsup > nctot)
    ncsup <- nctot
  
  if (center)
    centers <- colMeans(x)
  else
    centers <- rep(0,d)
  if (scale)
    scales <- apply(x,2,'sd')
  else
    scales <- rep(1,d)
  
  # this may produce NaNs in x if some columns have zero variance
  x <- t((t(x)-centers)/scales)
  
  # check if there are NA/NaNs in some columns, exclude those
  exclude <- c( exclude, which(is.na(colSums(x))) )
  smaller_than_mean <- rowSums( t(x) < colMeans(x) )
  exclude <- c(exclude, which(smaller_than_mean==1 | smaller_than_mean == n-1))
  ok <- setdiff(1:d, exclude)
  
  v_all <- matrix(0, nrow=d, ncol=nctot)
  b_all <- matrix(0, nrow=d, ncol=nctot)
  latent <- matrix(nrow=n, ncol=nctot)
  
  if (ncsup > 0) {
    
    print('Computing the supervised principal components..')
    
    for (k in 1:ncsup) {
      
      
      if (permtest) {
        if (permtest_type == 'max-marginal') {
          
          pval <- featscore.test(x,y,exclude=exclude,perms=perms,test.max=T)
          if (!any(pval < alpha)) {
            k <- k-1
            print(sprintf('SPC %d failed the max-marginal permutation test, so stopping supervised iteration.', k+1))
            break
          }
          
        } else if (permtest_type == 'marginal') {
          
          pval <- featscore.test(x,y,exclude=exclude,perms=perms)
          if (!any(pval < alpha)) {
            k <- k-1
            print(sprintf('SPC %d failed the marginal permutation test, so stopping supervised iteration.', k+1))
            break
          }
        }
      }
      
      # find the marginal score threshold so that the first PC maximizes the score
      if (is.factor(y)) {
        pcas <- lapply(levels(y), function(c) 
          spcs(x, y==c, thresh=thresh, nthresh=nthresh, exclude=exclude, window=window, nc=1, method=method) )
        scores <- vapply(pcas, function(pca) pca$score, 1.0)
        pca <- pcas[[which.max(scores)]]
      } else {
        pca <- spcs(x,y, thresh=thresh, nthresh=nthresh, exclude=exclude, window=window, nc=1, method=method)
      }
      
      
      if (permtest) {
        
        if (permtest_type == 'full') {
          fails <- 0
          for (rep in 1:perms) {
            scoreperm <- spcs(x, sample(y), thresh=thresh, nthresh=nthresh, 
                              exclude=exclude, window=window, nc=1, method=method)$score
            if (scoreperm > pca$score)
              fails <- fails+1
            if (fails/perms > alpha) {
              k <- k-1
              print(sprintf('SPC %d failed the permutation test, so stopping supervised iteration.', k+1))
              break
            }
          }
          if (fails/perms > alpha)
            break
        }
      }
      
      if (pca$score < min_score) {
        # score under the minimum, so terminate SPC computation
        k <- k-1
        print(sprintf('Could extract only %d supervised PCs although %s was asked.',k,ncsup))
        # warning('Could extract only %d supervised PCs although .')
        break
      }
      
      # the latent values and the principal component
      z <- as.vector(pca$x[,1])
      v <- as.vector(pca$rotation[,1])
      latent[,k] <- z
      
      # subract from x the variance explained by direction v
      b <- rep(0, d)
      b[ok] <- colSums(x[,ok,drop=F] * z) / sum(z^2)
      x <- x - t(t(matrix(z, nrow=n, ncol=d)) * b)
      
      # remove cols with very small variance for numerical stability
      # this usually occurs when v contains only one feature
      epsilon <- 1e-9
      xstd <- apply(x, 2, stats::sd)
      x[,xstd < epsilon] <- 0 
      
      v_all[,k] <- v
      b_all[,k] <- b
      
      # these are needed during the next iteration
      # exclude <- c(exclude, pca$ind) # should we exclude those variables already in some PC?
      # left <- setdiff(1:d,exclude)
      
      if (verbose) {
        print(sprintf('%d SPCs computed.', k, ncsup))
      }
    }
  } else {
    k <- 0
  }
  
  # ncsup may be smaller than initially specified if the iteration terminated
  # due to permutation test
  ncsup <- k
  
  if (!compute_unsup) {
    # nctot was initially unset indicating that we do not want to compute supervised components,
    # so remove the surplus from these arrays
    nctot <- ncsup
    v_all <- v_all[,1:nctot,drop=F]
    b_all <- b_all[,1:nctot,drop=F]
    latent <- latent[,1:nctot,drop=F]
  }
  
  if (ncsup < nctot) {
    # compute standard pc with the rest of the data variation
    print('Computing the unsupervised principal components..')
    if (length(ok) < nctot-ncsup)
      nctot <- length(ok)-ncsup # how many unsupervised PCs we can compute
    temp <- stats::prcomp(x[,ok], rank.=nctot-ncsup, center=F, scale.=F)
    v_all[ok,(ncsup+1):nctot] <- temp$rotation[, 1:(nctot-ncsup)]
    latent[,(ncsup+1):nctot] <- temp$x[, 1:(nctot-ncsup)]
  }
  
  # compute the rotation matrix from vs
  v <- v_all # store the original vs
  rotation <- v_all
  if (ncsup > 0 && nctot > 1) {
    for (j in seq(nctot, 2)) {
      for (h in seq(min(j-1,ncsup), 1))
        rotation[,j] <- rotation[,j] - sum(b_all[,h]*rotation[,j])*v_all[,h]
    }
  }
  
  if (normalize) {
    # scale the extracted features so that they have standard deviation of one
    sz <- apply(latent,2,'sd')
    latent <- t(t(latent)/sz)
    rotation <- t(t(rotation)/sz)
  }

  if (verbose)
    print('Done.')
  
  res <- list(w=rotation, z=latent, v=v, sdev=apply(latent,2,'sd'), 
              ncsup=ncsup, centers=centers, scales=scales, exclude=exclude)
  
  class(res) <- 'ispca'
  return(res)
  
}




#' @export
predict.ispca <- function(object, xnew, ...) {
  predict.dimred(object, xnew)
}

#' @export 
coeff.transform.ispca <- function(object, beta, alpha, ...) {
  coeff.transform.dimred(object, beta, alpha)
}


















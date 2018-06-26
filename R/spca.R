
#' Supervised (and unsupervised) principal components 
#'
#' Computes dimension reduction based on the supervised principal components
#' algorithm. In essense, algorithm performs a screening step based on univariate
#' scores for the features, and then computes standard PCA on using only the retained
#' features.
#' 
#' @param x The original feature matrix, columns denoting the features and rows the instances.
#' @param y A vector with the observed target values we try to predict using \code{x}.
#' Can be factor for classification problems. If missing, then this function computes
#' standard unsupervised principal components.
#' @param nctot Total number of latent features to extract.
#' @param ncsup Maximum number of latent features to extract that use supervision.
#' If \code{nctot > ncsup}, then the remaining \code{nctot-ncsup} features are computed in
#' unsupervised manner (or ignored if \code{sup.only=TRUE}).
#' @param window Maximum number of features that will survive the screening and from which
#' the supervised components are computed. Affects also how the \code{screenthresh}-argument
#' is interpreted.
#' @param exclude Columns (variables) in \code{x} to ignore when extrating the new features.
#' @param verbose Whether to print some messages along the way. 
#' @param normalize Whether to scale the extracted features so that they all have standard deviation
#'  of one.
#' @param preprocess Whether to center and scale the features before extracting the new features.
#' @param alpha Significance level for the p-values of the univariate scores used to determine
#'  which features survive the screening and are used to compute the supervised components.
#' @param perms Number of permutations to estimate the p-values for univariate scores.
#' @param screenthresh Value between 0 and 1 (or \code{NULL}). If not \code{NULL}, then no permutation tests are
#' run, and the supervised components
#' are computed among those features that have their univariate score equal or larger than this. 
#' Value 1 means that only the feature with the highest score survives the screening, whereas value 0 means that
#' the top \code{min(window, ncol(x))} survive the screening. Overwrites also \code{nfeat}-argument.
#' @param nfeat Number of features to retain in the screening step. If this option is used, then
#' the algorithm does not perform the permutation tests for the p-values, but instead computes the 
#' supervised components from those features that have their univariate score among the \code{nfeat} 
#' highest scores (in this case \code{perms} and \code{alpha} are ignored).
#' @param sup.only If \code{TRUE}, then no unsupervised components are ever computed even if the number of
#' supervised components that could be extracted was less than \code{nctot}.
#' @param ... Currently ignored.
#'
#'
#' @return spca-object that is similar to the object returned by \code{\link[stats]{prcomp}}.
#' The object will have the following elements:
#' \describe{
#'  \item{\code{w}}{The projection (or rotation) matrix W, that transforms the original data 
#'  \eqn{X} into the new features \eqn{Z = X W} .}
#'  \item{\code{z}}{The extracted latent features corresponding to the training inputs \eqn{X}.}
#'  \item{\code{v}}{Matrix \eqn{V} that is used to compute \eqn{W} when combining supervised and
#'  unsupervised components (see the Piironen and Vehtari (2018) for more information).}
#'  \item{\code{sdev}}{Standard deviations of the new features.}
#'  \item{\code{centers}}{Mean values for the original variables.}
#'  \item{\code{scales}}{Scales of the original variables.}
#'  \item{\code{exclude}}{Excluded variables.}
#' }
#' 
#' @section Details:
#' 
#' In the original paper, the authors proposed estimating the screening threshold
#' using cross-validation for the model obtained when the extracted features are used
#' for regression or classification. This implementation performs the screening
#' based on the estimated p-values for the univariate scores (these are estimated using
#' a permutation test) and the screening step retains only those features with p-value
#' less than the specified level \code{alpha}.
#' 
#' @section References:
#' 
#' Bair, E., Hastie, T., Paul, D., and Tibshirani, R. (2006).
#' Prediction by supervised principal components. \emph{Journal
#' of the American Statistical Association}, 101(473):119-137.
#' 
#' Piironen, J. and Vehtari, A. (2018). Iterative supervised principal components.
#' In \emph{Proceedings of the 21st International Conference on Artificial
#' Intelligence and Statistics (AISTATS) PMLR 84: 106-114}.
#'
#' @examples
#' \donttest{
#' ### 
#' dr <- spca(x,y, nctot=2)
#' z <- predict(dr, x) # the latent features
#' }
#'

#' @export
spca <- function(x, y=NULL, nctot=NULL, ncsup=NULL, window=500,
                 exclude=NULL, verbose=TRUE, normalize=FALSE,
                 preprocess=TRUE, alpha=NULL, perms=1000, 
                 screenthresh=NULL, nfeat=NULL, sup.only=FALSE, ...) {
  
  n <- NROW(x)
  d <- NCOL(x)
  
  if (!is.null(screenthresh)) {
    # a relative screening threshold given, so find out what is the number of features this
    # threshold corresponds to
    scores <- featscore(x,y,exclude=exclude)
    scores_sorted <- sort(scores, decreasing=T)
    smax <- scores_sorted[1]
    smin <- scores_sorted[min(length(scores),window)]
    nfeat <- sum(scores >= smin + screenthresh*(smax - smin)) 
  }
  
  if (is.null(alpha))
    alpha <- 0.001
  if (is.null(nctot)) {
    if (!is.null(nfeat))
      nctot <- min(nfeat, min(n-1,d))
    else 
      nctot <- min(n-1,d)
  }
  else if (nctot > min(n-1,d))
    stop('nctot cannot exceed min(n-1,d)')
  if (is.null(ncsup))
    ncsup <- nctot #ceiling(nctot/2)
  if (is.null(y)) {
    ncsup <- 0
  }
  
  if (preprocess) {
    # center and scale x for computation
    centers <- colMeans(x)
    scales <- apply(x,2,'sd')
    x <- scale(x) # this may produce NaNs in x if some columns have zero variance
  } else {
    centers <- rep(0,d)
    scales <- rep(1,d)
  }
  
  # check if there are NA/NaNs in some columns, exclude those
  exclude <- c( exclude, which(is.na(colSums(x))) )
  ok <- setdiff(1:d, exclude)
  
  v_all <- matrix(0, nrow=d,ncol=nctot)
  b_all <- matrix(0, nrow=d,ncol=nctot)
  latent <- matrix(nrow=n, ncol=nctot)
  rotation <- matrix(nrow=d, ncol=nctot)
  pval <- NULL
  
  
  # Procedure:
  # 1) compute the marginal scores
  # 2) threshold alpha for the permutation tests => compute ncsup SPCs
  # 3) subtract the variation explained by the SPCs from x
  # 4) compute nctot - ncsup unsupervised PCs
  # 5) compute the rotation matrices
  
  if (ncsup > 0) {
    
    if (verbose)
      print('Performing permutation tests for the marginal scores..')
    
    if (!is.null(nfeat)) {
      # screen out all except those among nfeat most relevant (based on the univariate scores)
      scores <- featscore(x,y, exclude=exclude)
      indsup <- order(scores, decreasing = T)[1:nfeat]
    } else {
      # permutation test for the marginal scores, find out subset for the supervised PCs
      pval <- featscore.test(x,y, exclude=exclude, perms=perms)
      indsup <- which(pval < alpha)
      if (length(indsup) > window) {
        # more than window significant features, so choose only the window most relevant ones
        scores <- featscore(x,y,exclude=exclude)
        indsup <- order(scores, decreasing = T)[1:window]
      }
    }
    
    if (length(indsup) > 0) {
      
      if (verbose)
        print('Computing the supervised principal components..')
      
      # supervised PCs
      spcs <- prcomp(x[, indsup], center=F, scale.=F)
      
      if ( NCOL(spcs$rotation) < ncsup )
        ncsup <- NCOL(spcs$rotation)
      
      v_all[indsup,1:ncsup] <- spcs$rotation[,1:ncsup]
      latent[,1:ncsup] <- spcs$x[,1:ncsup]
      
      # subtract the variation explained by the SPCs from x
      for (k in 1:ncsup) {
        
        z <- latent[,k]
        b <- rep(0, d)
        b[ok] <- colSums(x[,ok,drop=F] * z) / sum(z^2)
        x <- x - t(t(matrix(z, nrow=n, ncol=d)) * b)
        b_all[,k] <- b
      }
    } else {
      # none of the variables passed the permutation test
      ncsup <- 0
    }
  }
  
  rotation <- v_all
  
  if (!sup.only) {
    if (ncsup < nctot) {
      
      if (verbose)
        print('Computing the unsupervised principal components..')
      
      # compute standard pc with the rest of the data variation
      if (length(ok) < nctot-ncsup)
        nctot <- length(ok)-ncsup # how many unsupervised PCs we can compute
      pcs <- prcomp(x[,ok], center=F, scale.=F)
      v_all[ok,(ncsup+1):nctot] <- pcs$rotation[, 1:(nctot-ncsup)]
      latent[,(ncsup+1):nctot] <- pcs$x[, 1:(nctot-ncsup)]
      
      # compute the rotation matrix from vs
      # v <- v_all # store the original vs
      rotation <- v_all
      if (ncsup > 0 && nctot > 1) {
        for (j in seq(nctot, ncsup+1)) {
          for (h in seq(ncsup, 1))
            rotation[,j] <- rotation[,j] - sum(b_all[,h]*rotation[,j])*v_all[,h]
        }
      }
    }
  } else {
    # remove the NA-columns from the matrices that will be returned
    if (nctot > ncsup) {
      v_all <- v_all[,1:ncsup,drop=F]
      b_all <- b_all[,1:ncsup,drop=F]
      latent <- latent[,1:ncsup,drop=F]
      rotation <- rotation[,1:ncsup,drop=F]
    }
  }
  
  if (normalize) {
    sz <- apply(latent,2,'sd')
    ok <- sz > 1e-6 
    latent <- t(t(latent[,ok])/sz[ok])
    rotation <- t(t(rotation[,ok])/sz[ok])
  }
  
  if (verbose)
    print('Done.')
  
  res <- list(w=rotation, z=latent, v=v_all, sdev=apply(latent,2,'sd'), 
              pval=pval, centers=centers, scales=scales, exclude=exclude)
  
  class(res) <- 'spca'
  return(res)
  
}




#' @export
predict.spca <- function(object, xnew, ...) {
  predict.dimred(object, xnew)
}

#' @export 
coeff.transform.spca <- function(object, beta, alpha, ...) {
  coeff.transform.dimred(object, beta, alpha)
}











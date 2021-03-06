

#
# Function for computing the first PC using the power method.
# Might sometimes be faster than the SVD solution, but in general
# can have slow convergence if the first PC does not explain most of
# the variance in x.
#
power.pc <- function(x, maxiter = 10000, tol = 1e-4) {

  # find the first principal component of x using the power method
  d <- ncol(x)

  if (d == 1) {
    v <- 1
  } else {
    vprev <- stats::rnorm(d)
    for (i in 1:maxiter) {
      v <- t(x) %*% (x %*% vprev)
      v <- v / sqrt(sum(v * v))
      if (max(abs(v - vprev)) < tol) {
        break
      }
      vprev <- v
    }
    if (i == maxiter) {
      warning("Maximum number of iterations reached.")
    }
  }
  u <- x %*% v
  return(list(x = u, rotation = v))
}



#
# Auxiliary function to compute supervised PCs (screening + PCA)
#
spcs <- function(x, y, thresh = NULL, nthresh = NULL, exclude = NULL, nc = 1,
                 window = 100, optim_only = TRUE, method = "power", epsilon = 1e-6, ...) {
  # workhorse for computing supervised PCs
  # thresh is relative to the largest score statistic

  if (is.null(thresh)) {
    if (is.null(nthresh)) {
      nthresh <- 10
    }
    thresh <- seq(0, 1 - epsilon, len = nthresh)
  }

  if (any(thresh > 1 | thresh < 0)) {
    stop("thresh must be between 0 and 1")
  }

  # compute the scores for each feature that are used for screening (thresholding)
  D <- NCOL(x)
  scores <- featscore(x, y, exclude = exclude)
  scores[scores < epsilon] <- 0 # exclude very small scores for numerical stability
  max_score <- max(scores)
  cand <- order(scores, decreasing = T)
  cand <- cand[1:sum(scores > 0)] # only those x which have nonzero score
  if (length(cand) > window) {
    cand <- cand[1:window]
  }
  ncand <- length(cand)

  if (ncand == 0) {
    stop("Something went wrong: no variables with a nonzero score.")
  }

  if (length(thresh) > ncand) {
    # more thresholds than variables, so loop through all subset sizes
    subsets <- lapply(1:ncand, function(k) cand[1:k])
  } else {
    # those variables which have score above a certain threshold
    upper <- max(scores[cand])
    lower <- min(scores[cand])
    subsets <- lapply(
      (upper - lower) * thresh + lower,
      function(th) cand[scores[cand] >= th]
    )
  }

  pcas <- lapply(subsets, function(ind) {
    if (length(ind) == 0) {
      stop("Something went wrong: the active set became empty after screening.")
    } else {
      if (method == "power" && nc == 1) {
        # use the power method
        pca <- power.pc(x[, ind, drop = F])
      } else if (method == "svd") {
        # use svd method
        pca <- stats::prcomp(x[, ind, drop = F], ...)
        if (NCOL(pca$rotation) > nc) {
          pca$rotation <- pca$rotation[, 1:nc, drop = F]
          pca$x <- pca$x[, 1:nc, drop = F]
          pca$sdev <- pca$sdev[1:nc]
        }
      } else {
        stop(paste0("Got an unknown method '", method, "'."))
      }
      pca$ind <- ind
      v <- pca$rotation
      pca$rotation <- matrix(0, nrow = D, ncol = NCOL(v))
      pca$rotation[ind, ] <- v
      pca$score <- featscore(pca$x[, 1], y)
    }
    pca
  })

  if (optim_only) {
    scores <- vapply(pcas, function(pca) pca$score, 0.1)
    return(pcas[[which.max(scores)]])
  } else {
    if (length(thresh) == 1) {
      return(pcas[[1]])
    } else {
      return(pcas)
    }
  }
}




#' @export
predict.dimred <- function(object, xnew, ...) {
  # map xnew to the latent variable space znew
  d <- length(object$scales)
  if (is.vector(xnew)) {
    xnew <- matrix(xnew, ncol = d)
  }

  ok <- setdiff(1:d, object$exclude)
  xnew_standard <- t((t(xnew[, ok, drop = F]) - object$centers[ok]) / object$scales[ok])
  return(xnew_standard %*% object$w[ok, , drop = F])
}

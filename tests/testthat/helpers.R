
get.test.data <- function(n, sigma = 0.3, sigma_f = 0.1, latent_dim = 3, nfeat_per_latent = 2) {
  f <- matrix(rnorm(n * latent_dim), nrow = n)
  coef <- rep(0, latent_dim)
  coef[1:2] <- c(1, 0.5)
  func <- f %*% coef
  x <- matrix(nrow = n, ncol = nfeat_per_latent * latent_dim)
  xmat <- list()
  for (j in 1:latent_dim) {
    xmat[[j]] <- sigma_f * matrix(rnorm(n * nfeat_per_latent), nrow = n) + f[, j]
  }
  x <- do.call(cbind, xmat)
  y <- func + rnorm(n) * sigma
  return(list(x = x, y = y))
}



#' Transform linear model coefficients to original feature space
#'
#' Transforms linear model coefficients from a model fitted using latent features z
#' back to the original feature space x. This is possible since the dimension reduction methods
#' PCA, SPCA and ISPCA are all linear, and consequently, if a linear model is fitted using 
#' features z, it will correspond exactly to a linear model with features x.
#' This function can be used to compute the regression coefficients and intercept of that model.
#' @name coef_transform
#' 
#' @param object Dimension reduction object.
#' @param beta Regression coefficients in the z-space (latent space).
#' @param alpha Intercept in the z-space (latent space).
#' @param ... Currently ignored.
#'
#' @return A list with elements \code{beta} and \code{alpha} that give the model's 
#' regression coefficients and intercept, respectively, in the original feature space.
#' 
#'
#' @examples
#' \donttest{
#' ### 
#' 
#' # load data
#' data('ovarian', package = 'dimreduce')
#' x <- ovarian$x
#' y <- ovarian$y
#' 
#' # perform dimension reduction (use first two latent features)
#' dr <- ispca(x,y, nctot=2, normalize=TRUE)
#' z <- predict(dr,x)
#' 
#' # fit model
#' data <- data.frame(z, y)
#' if (requireNamespace('rstanarm', quietly=TRUE)) {
#' 
#'   # Bayesian logistic regression
#'   model <- rstanarm::stan_glm(y ~ ., data = data, family = binomial(), chains=2, iter=500)
#'   param_z <- as.data.frame(model)
#'   alpha_z <- t(param_z[,1])
#'   beta_z <- t(param_z[,2:3])
#'   
#'   # transform the linear model back to original feature space
#'   param_x <- coef_transform(dr, beta_z, alpha_z)
#'   beta_x <- param_x$beta
#'   alpha_x <- param_x$alpha
#'   
#'   intervals <- apply(beta_x, 1, quantile, probs=c(0.05,0.5,0.95))
#'   
#'   if (require('ggplot2', quietly = TRUE)) {
#'     # plot the coefficient posteriors in the original space
#'     p <- ncol(x)
#'     ggplot() + 
#'       geom_linerange(aes(x=1:p, ymin=intervals[1,], ymax=intervals[3,])) +
#'       geom_point(aes(x=1:p, y=intervals[2,])) +
#'       ylab('Coefficient') + 
#'       xlab('Feature')
#'   }
#'   
#' } else {
#' 
#'   # logistic regression
#'   model <- glm(y ~ ., data = data, family = binomial())
#'   param_z <- coef(model)
#'   alpha_z <- param_z[1]
#'   beta_z <- param_z[2:3]
#'   
#'   # transform the linear model back to original feature space
#'   param_x <- coef_transform(dr, beta_z, alpha_z)
#'   beta_x <- param_x$beta
#'   alpha_x <- param_x$alpha
#'   
#'   if (requireNamespace('ggplot2', quietly = TRUE)) {
#'     # plot the coefficients in the original space
#'     p <- ncol(x)
#'     ggplot2::ggplot() + 
#'       ggplot2::geom_point(ggplot2::aes(x=1:p, y=beta_x)) +
#'       ggplot2::ylab('Coefficient') + 
#'       ggplot2::xlab('Feature')
#'   }
#' }
#' 
#' }
#'
#'

#' @export
coef_transform.dimred <- function(object, beta, alpha, ...) {
  # transform linear regression coefficients from the z-space to 
  # the original x-space
  beta_x <- (object$w %*% beta)/object$scales
  alpha_x <- alpha - colSums(object$centers*beta_x)
  return(list(beta=beta_x, alpha=alpha_x))
}


#' @export
coef_transform <- function (object, ...) {
  UseMethod("coef_transform", object)
}





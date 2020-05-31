

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
#' data('prostate', package = 'dimreduce')
#' x <- prostate$x
#' y <- prostate$y
#'
#' # perform dimension reduction (use first two latent features)
#' dr <- ispca(x,y, nctot=2)
#' z <- predict(dr,x)
#'
#' # fit model
#' fit <- glm(y ~ z, family = binomial())
#' param_z <- coef(model)
#' alpha_z <- param_z[1]
#' beta_z <- param_z[2:3]
#'
#' # transform the linear model back to original feature space
#' param_x <- coef_transform(dr, beta, alpha)
#' beta_x <- param_x$beta
#' alpha_x <- param_x$alpha
#'
#' # plot the regression coefficients in the original space
#' plot(beta_x)
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





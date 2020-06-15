#' Dimension reduction for supervised learning
#'
#' @docType package
#' @name dimreduce
#'
#'
#' @description
#'
#' \pkg{dimreduce} is an R package that provided functions for
#' (supervised) dimension reduction.
#'
#' @section Functions:
#'
#' \describe{
#'  \item{\link{spca}, \link{ispca}}{
#'  Supervised PCA (SPCA) and iterative supervised PCA (ISPCA), that are useful techniques for
#'  dimension reduction. \link{spca} can also be used to  compute the standard unsupervised PCA.
#'  }
#'  \item{\link{coef_transform}}{
#'  Function that can be used to map linear model regression coefficients fitted using
#'  latent features from SPCA or ISPCA back to the original feature space.
#'  This can be useful for analysing the model.
#'  }
#'  \item{\link{featscore}, \link{featscore.test}}{
#'  Functions for computing univariate relevance scores and significance tests for the features,
#'  which can be used for screening.
#'  }
#' }
#'
#' @section References:
#'
#' Bair, E., Hastie, T., Paul, D., and Tibshirani, R. (2006).
#' Prediction by supervised principal components. \emph{Journal
#' of the American Statistical Association}, 101(473):119-137.
#'
#' Neal, R. and Zhang, J. (2006). High dimensional classification with Bayesian
#' neural networks and Dirichlet diffusion trees.
#' In Guyon, I., Gunn, S., Nikravesh, M., and Zadeh, L. A., editors,
#' \emph{Feature Extraction, Foundations and Applications}, pages 265-296. Springer.
#'
#' Piironen, Juho and Vehtari, Aki (2018). Iterative supervised principal components.
#' In \emph{Proceedings of the 21st International Conference on Artificial
#' Intelligence and Statistics (AISTATS) PMLR 84: 106-114}.
#'
NULL

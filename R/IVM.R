#' Perform IVM selection for sparse GP regression
#'
#' @param predictors a numeric matrix of predictors (X)
#' @param activeSetSize the number of points to select for the active set
#' @param variance the variance of the Gaussian noise function
#' @param kernel.function a positive definite function which takes two vectors
#'   x_1 and x_2 and returns the covariance between the two points
#'
#' @return An object of class IVM. Active set indices are stored in activeSet
#' @export
IVM.regression <- function(predictors,
                           activeSetSize,
                           variance,
                           kernel.function) {
  predictors <- as.matrix(predictors)
  numSamples <- nrow(predictors)
  out <- list()
  class(out) <- "IVM"
  activeSet <- c()
  inactiveSet <- seq(numSamples)
  M <- matrix(0, activeSetSize, numSamples)
  cov.diag <- numeric(numSamples)
  for (i in seq_along(cov.diag)) {
    cov.diag[i] <- kernel.function(predictors[i,], predictors[i, ])
  }

  for (i in 1:activeSetSize) {
    max.j.index <- NA
    max.entropyDiff <- -Inf
    max.nu <- NA
    for (j.index in seq_along(inactiveSet)) {
      j <- inactiveSet[j.index]
      nu <- 1 / (cov.diag[j] + variance)
      entropyDiff <- -1 / 2 * log(1 - nu * cov.diag[j])
      if (entropyDiff > max.entropyDiff) {
        max.j.index <- j.index
        max.entropyDiff <- entropyDiff
        max.nu <- nu
      }
    }
    max.j <- inactiveSet[max.j.index]

    K_j <- numeric(numSamples)
    for (k in 1:numSamples) {
      K_j[k] <- kernel.function(predictors[max.j, ], predictors[k, ])
    }
    s_j <- K_j - t(M) %*% M[, max.j]

    M[i, ] <- sqrt(max.nu) * s_j

    cov.diag <- cov.diag - max.nu * diag(s_j %*% t(s_j))

    activeSet <- c(activeSet, max.j)
    inactiveSet <- inactiveSet[-max.j.index]

  }
  out$activeSet <- activeSet
  out$M <- M
  out$cov.diag <- cov.diag
  return(out)
}

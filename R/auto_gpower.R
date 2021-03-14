source("R/gpower.R")
#' Find GPower rho for proportion of sparcity
#'
#' Using binary search find the value of rho for \code{\link{gpower}} for which
#' the the proporting of values equal to zero is closest to the value of
#' prop_sparse. It forwards the other settings to the \code{\link{gpower}}
#' function.
#'
#' @param prop_sparse the percentage of the total values of the loadings matrix
#' which is equal to zero.
#' @inheritParams gpower
#'
#' @export
auto_gpower <- function(X,
                        k,
                        prop_sparse,
                        penalty = c("l0", "l1"),
                        block = c(TRUE, FALSE),
                        mu = NA,
                        iter_max = 1000,
                        epsilon = 1e-4) {
  UseMethod("auto_gpower")
}


#' @export
auto_gpower.default <- function(X,
                                k,
                                prop_sparse,
                                penalty = "l1",
                                block = FALSE,
                                mu = NA,
                                iter_max = 1000,
                                epsilon = 1e-4) {
  # Tunes rho to get desired proportion of sparsity

  n <- ncol(X)
  n_zeros_X <- floor(n * k * prop_sparse)
  n_zeros_sparse <- 0

  if (n_zeros_X == 0) {
    # No sparsity
    gpower(X, k, 0, penalty, block, mu, iter_max, epsilon)
  }

  # Starting bounds of binary search
  lower <- 0
  if (penalty == "l0") {
    upper <- max(norm(X, "2"))
  }
  if (penalty == "l1") {
    upper <- max(norm(X, "2"))^2
  }
  i <- 0

  while (n_zeros_X != n_zeros_sparse & iter_max > i) {
    cut <- (lower + upper) / 2

    if (penalty == "l0") {
      gamma <- (cut^2)
    }
    if (penalty == "l1") {
      gamma <- cut
    }

    Z <- gpower(X, k, gamma, penalty, block, mu, iter_max, epsilon)

    n_zeros_sparse <- sum(rowSums(Z$loadings == 0))

    if (n_zeros_X > n_zeros_sparse) {
      lower <- gamma
    }
    if (n_zeros_X < n_zeros_sparse) {
      upper <- gamma
    }

    i <- i + 1
  }
  cat(
    "After",
    i,
    "iterations, rho",
    gamma,
    "achieves",
    prop_sparse,
    "proportion of sparseness",
    sep = " "
  )

  Z
}


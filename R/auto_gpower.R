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
#' @examples
#' set.seed(360)
#' p <- 20
#' n <- 50
#' k <- 5
#' data <- matrix(stats::rnorm(p * n), nrow = p, ncol = n)
#' prop_sparse <- 0.1
#' mu <- 1
#'
#' auto_gpower(data, k, prop_sparse, penalty = 'l1', center = TRUE, block = FALSE)
#'
#' @export
auto_gpower <- function(data, k, prop_sparse, penalty = c("l0", "l1"), center = c(TRUE, FALSE), block = c(TRUE, FALSE), mu = 1, iter_max = 1000, epsilon = 1e-04) {
    UseMethod("auto_gpower")
}


#' @export
auto_gpower.default <- function(data, k, prop_sparse, penalty = "l1", center = TRUE, block = FALSE, mu = 1, iter_max = 1000, epsilon = 1e-04) {
    # Tunes rho to get desired proportion of sparsity

    n <- ncol(data)
    n_zeros_data <- floor(n * k * prop_sparse)
    n_zeros_sparse <- 0

    if (n_zeros_data == 0) {
        # No sparsity
        gpower(data, k, 0, penalty, center, block, mu, iter_max, epsilon)
    }

    # Starting bounds of binary search
    lower <- 0
    if (penalty == "l0") {
        upper <- max(norm(data, "2"))
    }
    if (penalty == "l1") {
        upper <- max(norm(data, "2"))^2
    }
    i <- 0

    while (n_zeros_data != n_zeros_sparse & iter_max > i) {
        cut <- (lower + upper)/2

        if (penalty == "l0") {
            gamma <- (cut^2)
        }
        if (penalty == "l1") {
            gamma <- cut
        }
        Z <- tryCatch(gpower(data = data, k = k, rho = x, penalty = penalty, center = center, block = block, mu = mu), error = function(e) {
            NULL
        })

        if (is.list(Z)) {
            n_zeros_sparse <- sum(rowSums(Z$loadings == 0))

            if (n_zeros_data > n_zeros_sparse) {
                lower <- gamma
            }
            if (n_zeros_data < n_zeros_sparse) {
                upper <- gamma
            }
        } else {
            upper <- gamma
        }



        i <- i + 1
    }
    cat("After", i, "iterations, rho", gamma, "achieves", prop_sparse, "proportion of sparseness", sep = " ")

    Z
}

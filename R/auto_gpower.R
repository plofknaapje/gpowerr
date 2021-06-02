#' Find GPower rho for proportion of sparcity
#'
#' Using binary search find the value of rho for \code{\link{gpower}} for which
#' the the proporting of values equal to zero is closest to the value of
#' prop_sparse. It forwards the other settings to the \code{\link{gpower}}
#' function.
#'
#' @param prop_sparse The percentage of the total values of the weights matrix
#' which is equal to zero.
#' @param accuracy The amount of digits to which to round prop_sparse and the
#'   sparsity of the weights.
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
auto_gpower <-
    function(data,
             k,
             prop_sparse,
             accuracy = 2,
             penalty = c("l0", "l1"),
             center = c(TRUE, FALSE),
             block = c(TRUE, FALSE),
             mu = 1,
             iter_max = 1000,
             epsilon = 1e-04) {
        UseMethod("auto_gpower")
    }


#' @export
auto_gpower.default <-
    function(data,
             k,
             prop_sparse,
             accuracy = 2,
             penalty = "l1",
             center = TRUE,
             block = FALSE,
             mu = 1,
             iter_max = 1000,
             epsilon = 1e-04) {
        # Tunes rho to get desired proportion of sparsity

        n <- ncol(data)
        prop_zeros_needed <- round(prop_sparse, accuracy)
        prop_zeros_current <- 0

        if (prop_zeros_needed == 0) {
            # No sparsity
            gpower(data, k, 0, penalty, center, block, mu, iter_max, epsilon)
        }

        if (block) {
            max_iterations <- 500
        } else {
            max_iterations <- 100
        }

        # Starting bounds of binary search
        lower <- 0
        upper <- 1
        i <- 0

        while (prop_zeros_needed != prop_zeros_current &
               max_iterations > i) {
            middle <- (lower + upper) / 2

            Z <- tryCatch(
                gpower(
                    data = data,
                    k = k,
                    rho = middle,
                    penalty = penalty,
                    center = center,
                    block = block,
                    mu = mu,
                    iter_max = iter_max,
                    epsilon = epsilon
                ),
                error = function(e) {
                    NULL
                }
            )

            if (is.list(Z)) {
                prop_zeros_current <-
                    round(sum(rowSums(Z$weights == 0)) / (n * k), accuracy)

                if (prop_zeros_needed > prop_zeros_current) {
                    lower <- middle
                }
                if (prop_zeros_needed < prop_zeros_current) {
                    upper <- middle
                }
                if (prop_zeros_needed == prop_zeros_current) {
                    break
                }
            } else {
                upper <- middle
            }
            i <- i + 1
        }

        if (is.list(Z)) {
            cat(
                "After",
                i,
                "iterations, rho",
                round(middle, 5),
                "achieves",
                round(sum(rowSums(
                    Z$weights == 0
                )) / (n * k), accuracy),
                "Proportion of Sparseness\n\n",
                sep = " "
            )

            Z
        } else {
            warning("Unable to find a suitable rho for specified prop_sparse")
        }
    }

#' Compute Sparse PCA using GPower method
#'
#' Generalized power method for sparse principal component analysis. Implements
#' the method developed by Journee et al. (2010) with a choice between a L1 and
#' L0 regularisation and a column based and block approach.
#'
#' @description{ \loadmathjax GPower uses four different optimization procedures for the four
#'   combinations between  \mjseqn{l_0} and \mjseqn{l_1} regularisation and single-unit or block
#'   computation. The function tries to find a weights matrix
#'   \mjseqn{W \in R^{n \times k}} which has the highest possible explained
#'   variance of the data matrix \mjseqn{X \in R^{p \times n}} under the regularisation
#'   constraints of the case. The matrix \mjseqn{Z \in R^{n \times k}} is used by some
#'   of the methods as an intermediate solution. Lambda is calculated by
#'   multiplying rho with the maximum possible value of lambda.
#'
#'   The objective function of the single unit case with \mjseqn{l_1} regularisation is
#'   \mjsdeqn{\hat{w} = \underset{\| w \| = 1}{\textrm{argmax}}{\| Xw \| - \lambda \| w \|_1}}
#'   For the single-unit case with the \mjseqn{l_0} regularisation, the objective function is
#'   \mjsdeqn{\hat{w} = \underset{\| z \| = 1}{\textrm{argmax}}\;
#'   \underset{\| w \| = 1}{\textrm{argmax}}{ (z^\top X w)^2 - \lambda \| w \|_0},}
#'   where the results are squared before gamma is subtracted instead of after.
#'   In order to compute more than 1 component, the matrix \mjseqn{X} is adjusted after each new component.
#'
#'   For the block cases, the following functions are used. For the case with \mjseqn{l_1} regularisation,
#'   \mjsdeqn{\hat{W} = \underset{Z \in M^p_k}{\textrm{argmax}} \sum_{j=1}^k
#'   {\underset{\| W_j \| = 1}{\textrm{argmax}} \mu_j Z_j^{\top} XW_j -
#'   \lambda_j {\| Z_j \|} }}
#'   and for the \mjseqn{l_0} regularisation case,
#'   \mjsdeqn{\hat{W} = \underset{Z \in M^p_k}{\textrm{argmax}} \sum_{j=1}^k
#'   \underset{\| \hat{W}_j \| = 1}{\textrm{argmax}} (\mu_j Z_i^{\top}
#'   X W_j)^2 - \lambda_j\| W_j \|_0}
#'
#'   All of these functions are optimized using the generalized power approach
#'   as described in the paper by Journée et al. (2010).}
#'
#' @param data Input matrix of size (p x n) with p < n.
#' @param k Number of components, 0 < k < p.
#' @param rho Relative sparsity weight factor of the optimization. Either a
#'   vector of floats of size k or float which will be repeated k times. 0 < rho
#'   < 1.
#' @param reg regularisation type to use in the optimization. Either 'l0' or 'l1'.
#'   The default is 'l1' since it performed best in experiments.
#' @param center Centers the data. Either TRUE or FALSE. The default is TRUE.
#' @param block Optimization method. If FALSE, the components are calculated
#'   individually. If TRUE, all components are calculated at the same time.
#'   The default is FALSE.
#' @param mu Mean to be applied to each component in the block. Either a vector
#'   of float of size k or a float which will be repeated k times. Only used if
#'   block is TRUE. The default is 1.
#' @param iter_max Maximum iterations when adjusting components with gradient
#'   descent. The default is 1000.
#' @param epsilon Epsilon of the gradient descent stopping function.
#'   The default is 1e-4.
#'
#' @return List containing: \describe{ \item{weights}{The PCA components}
#'   \item{scores}{Scores of the components on data} \item{a_approx}{Reconstructed
#'   version of data using the components} \item{prop_sparse}{Proportion of
#'   sparsity of the components} \item{exp_var}{Explained ratio of variance of
#'   the components} \item{centers}{Centers of matrix data if center == TRUE} }
#'
#' @examples
#' set.seed(360)
#' p <- 20
#' n <- 50
#' k <- 5
#' data <- matrix(stats::rnorm(p * n), nrow = p, ncol = n)
#' rho <- 0.1
#' # rho <- c(0.1, 0.2, 0.1, 0.2, 0.1)
#' mu <- 1
#' # mu <- c(1, 1.5, 0.5, 2, 1)
#'
#' # Single unit with l1 regularisation
#' gpower(data, k, rho, 'l1', TRUE)
#'
#' # Single unit with l0 regularisation
#' gpower(data, k, rho, 'l0', TRUE)
#'
#' # Block with l1 regularisation
#' gpower(data, k, rho, 'l1', TRUE, TRUE, mu)
#'
#' # Block with l0 regularisation
#' gpower(data, k, rho, 'l0', TRUE, TRUE, mu)
#'
#' @references Journee, M., Nesterov, Y., Richtarik, P. and Sepulchre, R. (2010)
#' Generalized Power Method for Sparse Principal Component Analysis. *Journal of
#' Machine Learning Research. 11*, 517-553.
#'
#' @export
gpower <-
    function(data,
             k,
             rho,
             reg = c("l0", "l1"),
             center = c(TRUE, FALSE),
             block = c(TRUE, FALSE),
             mu = 1 ,
             iter_max = 1000,
             epsilon = 1e-04) {
        UseMethod("gpower")
    }

#' @export
gpower.default <-
    function(data,
             k,
             rho,
             reg = "l1",
             center = TRUE,
             block = FALSE,
             mu = 1,
             iter_max = 1000,
             epsilon = 1e-04) {
        data <- as.matrix(data)

        sparsity_error <-
            "Sparcity is set too high, all entries of loading vector are zero"
        reg_method <- "Regularisation method not recognized"
        iter_warning <- "Maximum number of iterations was reached"

        # Checks ------------------------------------------------------------------

        if (any(is.na(data))) {
            warning("Missing values are omitted: na.omit(X).")
            X <- stats::na.omit(X)
        }

        if (any(is.na(rho))) {
            warning("NA found in rho")
        }

        if (block == 1 & any(is.na(mu))) {
            warning("Mu needs to be defined for block algoritm")
        }

        picked <- FALSE

        # Initialize gpower object ------------------------------------------------
        gpowerObj <- list(weights = NULL, scores = NULL)

        p <- nrow(data)
        n <- ncol(data)

        iter_max <- iter_max
        epsilon <- epsilon

        Z <- matrix(0, n, k, dimnames = list(colnames(data)))

        W <- matrix(0, n, k, dimnames = list(colnames(data)))

        if (length(rho) == 1) {
            rho <- rep(rho, k)
        }

        if (length(mu) == 1 & !any(is.na(mu))) {
            mu <- rep(mu, k)
        }
        # Add centering

        if (center) {
            data <- scale(data, scale = FALSE)
            gpowerObj$centers <- attr(data, "scaled:center")
        }

        X <- data

        # Single unit algorithm ---------------------------------------------------

        if (!picked & (k == 1 | (k > 1 & !block))) {
            if (reg == "l1") {
                for (c in 1:k) {
                    # Loop on the components
                    norm_a_i <- rep(0, n)

                    for (i in 1:n) {
                        norm_a_i[i] <- norm(X[, i], "2")
                    }

                    rho_max <- max(norm_a_i)
                    i_max <- which.max(norm_a_i)
                    rho_c <- rho[c] * rho_max

                    # Initialisation Allow for warm start
                    x <- X[, i_max] / norm_a_i[i_max]
                    f <- rep(0, iter_max)
                    iter <- 1

                    while (TRUE) {
                        X_x <- t(X) %*% x
                        t_resh <-
                            sign(X_x) * pmax(abs(X_x) - rho_c, 0)

                        # Cost function
                        f[iter] <- sum(t_resh ^ 2)

                        if (f[iter] == 0) {
                            # Sparsity is too high
                            warning(sparsity_error)
                            break
                        } else {
                            gradient <- X %*% t_resh
                            x <- gradient / norm(gradient, "2")
                        }

                        if (iter > 2) {
                            # Stopping criteria
                            if ((f[iter] - f[iter - 1]) / f[iter - 1] < epsilon |
                                iter > iter_max) {
                                if (iter > iter_max) {
                                    print(iter_warning)
                                    break
                                }
                                break
                            }
                        }

                        iter <- iter + 1
                    }

                    X_x <- t(X) %*% x
                    pattern <-
                        (abs(X_x) - rho_c) > 0  # Pattern of sparsity
                    z <- sign(X_x) * max(abs(X_x) - rho_c, 0)

                    if (max(abs(z > 0))) {
                        z <- z / norm(z, "2")
                    }

                    z <- pattern_filling(X, pattern, z)

                    # Deflate
                    y <- X %*% z
                    X <- X - y %*% z

                    W[, c] <- z
                }
            } else if (reg == "l0") {
                for (c in 1:k) {
                    # Loop on the components
                    rho_c <- rho[c]
                    norm_a_i <- rep(0, n)

                    for (i in 1:n) {
                        norm_a_i[i] <- norm(X[, i], "2")
                    }

                    rho_max <- max(norm_a_i)
                    i_max <- which.max(norm_a_i)
                    rho_c <- rho_c * rho_max ^ 2

                    # Initialisation
                    x <- X[, i_max] / norm_a_i[i_max]
                    f <- rep(0, iter_max)
                    iter <- 1

                    while (TRUE) {
                        X_x <- t(X) %*% x
                        t_resh <- pmax(X_x ^ 2 - rho_c, 0)

                        # Cost function
                        f[iter] <- sum(t_resh)

                        if (f[iter] == 0) {
                            # Sparsity is too high
                            warning(sparsity_error)
                            break
                        } else {
                            gradient <- X %*% ((t_resh > 0) * X_x)
                            x <- gradient / norm(gradient, "2")
                        }

                        if (iter > 2) {
                            # Stopping criteria
                            if ((f[iter] - f[iter - 1]) / f[iter - 1] < epsilon |
                                iter > iter_max) {
                                if (iter > iter_max) {
                                    print(iter_warning)
                                    break
                                }
                                break
                            }
                        }

                        iter <- iter + 1
                    }

                    pattern <-
                        ((t(X) %*% x) ^ 2 - rho_c > 0)  # Pattern of sparsity
                    y <- x
                    pattern_inv <- pattern == 0

                    z <- t(X) %*% y
                    z[pattern_inv] <- 0
                    norm_z <- norm(z, "2")
                    z <- z / norm_z
                    y <- y * norm_z

                    # Deflate
                    X <- X - y %*% t(z)

                    W[, c] <- z
                }
            } else {
                warning(reg_method)
            }
            picked <- TRUE
        }

        # Block algorithm with mu = 1 ---------------------------------------------
        if (!picked & (k > 1 & block & sum(mu == 1) == length(mu))) {
            norm_a_i <- rep(0, n)

            for (i in 1:n) {
                norm_a_i[i] <- norm(X[, i], "2")
            }

            i_max <- which.max(norm_a_i)

            # Initialization

            qr_decomp <-
                qr(cbind(X[, i_max] / norm_a_i[i_max], matrix(stats::rnorm(p * (
                    k - 1
                )), nrow = p)), LAPACK = T)  # LAPACK to get MatLab qr(X, 0) results

            x <- qr.Q(qr_decomp)
            rho_max <- qr.R(qr_decomp)


            f <- rep(0, iter_max)
            iter <- 1

            if (reg == "l1") {
                rho <- rho %*% rho_max

                while (TRUE) {
                    X_x <- t(X) %*% x
                    t_resh <-
                        pmax(abs(X_x) - kronecker(matrix(1, n, 1), rho), 0)
                    f[iter] <- sum(t_resh ^ 2)

                    if (f[iter] == 0) {
                        warning(sparsity_error)
                        break
                    } else {
                        gradient <- matrix(0, p, k)

                        for (i in 1:k) {
                            pattern <- t(t_resh[, i]) > 0
                            gradient[, i] <- tryCatch(
                                X[, pattern] %*% (t_resh[pattern, i] * sign(X_x[pattern, i])),
                                error = function(e) {
                                    NA
                                })
                            if (any(is.na(gradient[, i]))) {
                                stop("Value of rho is too high")
                            }

                        }

                        svd_decomp <- svd(gradient)
                        x <- svd_decomp$u %*% t(svd_decomp$v)
                    }

                    if (iter > 2) {
                        if (iter > iter_max) {
                            print(iter_warning)
                            break
                        }
                        if ((f[iter] - f[iter - 1]) / f[iter - 1] < epsilon) {
                            break
                        }
                    }
                    iter <- iter + 1
                }

                X_x <- t(X) %*% x

                for (i in 1:k) {
                    Z[, i] <- sign(X_x[, i]) * max(abs(X_x[, i]) - rho[i], 0)
                    if (max(abs(Z[, i]) > 0) > 0) {
                        Z[, i] <- Z[, i] / norm(Z[, i], "2")
                    }
                }

                pattern <-
                    abs(X_x) - kronecker(matrix(1, n, 1), rho) > 0
                W <- pattern_filling(X, pattern, Z, mu)

            } else if (reg == "l0") {
                rho <- rho %*% (rho_max ^ 2)

                while (TRUE) {
                    X_x <- t(X) %*% x
                    t_resh <-
                        pmax(X_x ^ 2 - kronecker(matrix(1, n, 1), rho), 0)
                    f[iter] <- sum(sum(t_resh))
                    if (f[iter] == 0) {
                        warning(sparsity_error)
                        break
                    } else {
                        gradient <- matrix(0, p, k)

                        for (i in 1:k) {
                            pattern <- t(t_resh[, i]) > 0

                            gradient[, i] <- tryCatch(
                                X[, pattern] %*% X_x[pattern, i],
                                error = function(e) {
                                    NA
                                })
                            if (any(is.na(gradient[, i]))) {
                                stop("Value of rho is too high")
                            }


                        }

                        svd_decomp <- svd(gradient)
                        x <- svd_decomp$u %*% t(svd_decomp$v)
                    }

                    if (iter > 2) {
                        if (iter > iter_max) {
                            print(iter_warning)
                            break
                        }
                        if ((f[iter] - f[iter - 1]) / f[iter - 1] < epsilon) {
                            break
                        }
                    }
                    iter <- iter + 1
                }

                pattern <-
                    (X_x ^ 2 - kronecker(matrix(1, n, 1), rho)) > 0
                pattern_inv <- pattern == 0
                W <- X_x
                W[pattern_inv] <- 0
                norm_z <- rep(0, k)

                for (i in 1:k) {
                    norm_z[i] <- norm(W[, i], "2")
                    if (norm_z[i] > 0) {
                        W[, i] <- W[, i] / norm_z[i]
                    }
                }
            } else {
                warning(reg_method)
            }
            picked <- TRUE
        }

        # Block algorithm with mu != 1 ---------------------------------------------
        if (!picked & (k > 1 & block & sum(mu == 1) < length(mu))) {
            norm_a_i <- rep(0, n)

            for (i in 1:n) {
                norm_a_i[i] <- norm(X[, i], "2")
            }

            i_max <- which.max(norm_a_i)

            # Initialization

            qr_decomp <-
                qr(cbind(X[, i_max] / norm_a_i[i_max], matrix(stats::rnorm(p * (
                    k - 1
                )), nrow = p)), LAPACK = T)  # LAPACK to get MatLab qr(X, 0) results

            x <- qr.Q(qr_decomp)
            rho_max <- qr.R(qr_decomp)


            f <- rep(0, iter_max)
            iter <- 1

            if (reg == "l1") {
                rho <- rho * mu %*% rho_max

                while (TRUE) {
                    X_x <- t(X) %*% x

                    for (i in 1:k) {
                        X_x[, i] <- X_x[, i] * mu[i]
                    }

                    t_resh <-
                        pmax(abs(X_x) - kronecker(matrix(1, n, 1), rho), 0)
                    f[iter] <- sum(t_resh ^ 2)

                    if (f[iter] == 0) {
                        warning(sparsity_error)
                        break
                    } else {
                        gradient <- matrix(0, p, k)

                        for (i in 1:k) {
                            pattern <- t(t_resh[, i]) > 0
                            gradient[, i] <-
                                X[, pattern] %*% (t_resh[pattern, i] * sign(X_x[pattern, i])) * mu[i]
                        }

                        svd_decomp <- svd(gradient)
                        x <- svd_decomp$u %*% t(svd_decomp$v)
                    }

                    if (iter > 2) {
                        if (iter > iter_max) {
                            print(iter_warning)
                            break
                        }
                        if ((f[iter] - f[iter - 1]) / f[iter - 1] < epsilon) {
                            break
                        }
                    }
                    iter <- iter + 1
                }

                X_x <- t(X) %*% x

                for (i in 1:k) {
                    X_x[, i] <- X_x[, i] * mu[i]
                    Z[, i] <-
                        sign(X_x[, i]) * max(abs(X_x[, i]) - rho[i], 0)
                    if (max(abs(Z[, i]) > 0) > 0) {
                        Z[, i] <- Z[, i] / norm(Z[, i], "2")
                    }
                }

                pattern <-
                    abs(X_x) - kronecker(matrix(1, n, 1), rho) > 0
                W <- pattern_filling(X, pattern, Z, mu)
            } else if (reg == "l0") {
                rho <- rho * (mu %*% rho_max) ^ 2

                while (TRUE) {
                    X_x <- t(X) %*% x

                    for (i in 1:k) {
                        X_x[, i] <- X_x[, i] * mu[i]
                    }

                    t_resh <-
                        pmax(X_x ^ 2 - kronecker(matrix(1, n, 1), rho), 0)
                    f[iter] <- sum(sum(t_resh))
                    if (f[iter] == 0) {
                        warning(sparsity_error)
                        break
                    } else {
                        gradient <- matrix(0, p, k)

                        for (i in 1:k) {
                            pattern <- t(t_resh[, i]) > 0
                            gradient[, i] <-
                                X[, pattern] %*% X_x[pattern, i]
                        }

                        svd_decomp <- svd(gradient)
                        x <- svd_decomp$u %*% t(svd_decomp$v)
                    }

                    if (iter > 2) {
                        if (iter > iter_max) {
                            print(iter_warning)
                            break
                        }
                        if ((f[iter] - f[iter - 1]) / f[iter - 1] < epsilon) {
                            break
                        }
                    }
                    iter <- iter + 1
                }

                X_x <- t(X) %*% x
                for (i in 1:k) {
                    X_x[, i] <- X_x[, i] * mu[i]
                }

                pattern <-
                    (X_x ^ 2 - kronecker(matrix(1, n, 1), rho)) > 0
                pattern_inv <- pattern == 0
                W <- t(X) %*% x
                W[pattern_inv] <- 0
                norm_z <- rep(0, k)

                for (i in 1:k) {
                    norm_z[i] <- norm(W[, i], "2")
                    if (norm_z[i] > 0) {
                        W[, i] <- W[, i] / norm_z[i]
                    }
                }
            } else {
                warning(reg_method)
            }
            picked <- TRUE
        }


        # Update gpower object ----------------------------------------------------
        scores <- data %*% W
        P <- t(data) %*% MASS::ginv(t(scores))
        data_approx <- scores %*% t(P)
        colnames(W) <- as.character(1:k)

        gpowerObj$weights <- W
        gpowerObj$scores <- scores
        gpowerObj$a_approx <- data_approx
        gpowerObj$prop_sparse <- sum(rowSums(W == 0)) / (n * k)

        # Add component variance (p22)
        AdjVar <- qr.R(qr(scores)) ^ 2

        ### Extract values on the diagonal
        comp_var <- rep(NA, k)
        for (i in 1:k) {
            comp_var[i] <- AdjVar[i, i]
        }

        # Explained variance ratio
        gpowerObj$comp_var <- sort(comp_var, decreasing = TRUE)
        gpowerObj$exp_var <-
            1 - (norm((data - data_approx), type = "F") / norm(data, type = "F")) ^
            2

        class(gpowerObj) <- "gpower"

        gpowerObj
    }


#' Prints the proportion of explained variance and the weights of the gpower
#' object.
#' @param x gpower object
#' @param print_zero_rows If FALSE, then only the rows which contain at least 1
#' non-zero value will be printed.
#' @param ... print() parameters
#' @export
print.gpower <- function(x, print_zero_rows = TRUE, ...) {
    # Print gpower --------------------------------------------------------------

    cat("Proportion of Explained Variance\n")
    print(round(x$exp_var, 3))
    cat("\nSparse weights:\n")
    if (print_zero_rows) {
        print(round(x$weights, 3))
    }
    else {
        row_has_nonzero <- apply(x$weights, 1, function(x){any(x != 0)})
        print(round(x$weights[row_has_nonzero, ], 3))
    }

}

#' Prints the proportion of explained variance and the proportion of sparseness.
#' @param object gpower object
#' @param ... summary() parameters
#' @export
summary.gpower <- function(object, ...) {
    # Summary gpower ------------------------------------------------------------
    x <-
        t(data.frame(
            var = round(object$exp_var, 3),
            prop_sparse = round(object$prop_sparse, 3)
        ))

    rownames(x) <-
        c("Proportion of Explained Variance", "Proportion of Sparsity")
    x
}

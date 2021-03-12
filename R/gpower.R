source("R/pattern_filling.R")

#' Compute Sparse PCA using GPower method
#'
#' Generalized power method for sparse principal component analysis. Implements
#' the method developed by Journee et al. (2010) with a choice between a L1 and
#' L0 penalty and a column based and block approach.
#'
#' GPower uses four different optimization procedures for the four combinations
#' between l0 and l1 penalty and single-unit or block computation. With the l0
#' penalty, the cardinality of the solutions is penalized. The objective
#' function of the single unit case with l1 penalty is
#' \denq{\phi_{l_{1}}^{2}(\gamma) = \max_{x\in S^{p}}\sum_{i=1}^{n}{\left [
#' \left |  a_{i}^{T}x\right | - \gamma  \right ]_{+}^{2}}}
#' where x represents the components and gamma is the penalty.
#' For the single unit case with the l0 penalty, the following function is used
#' \denq{\phi_{l_{0}}^{2}(\gamma) = \max_{x\in S^{p}}\sum_{i=1}^{n}{\left [
#' \left ( a_{i}^{T} x \right )^{2} - \gamma  \right ]_{+}}}
#' Where the results are squared before gamma is subtracted instead of after.
#' For the block cases, the following functions are used. First the case with l1
#' penalty
#' \denq{\phi_{l_{1},m}^{2}(\gamma) = \max_{x\in S^{p}}\sum_{j=1}^{m}{\sum_{i=1}
#' ^{n}{\left [ \mu_{j} \left | a_{i}^{T}x_{j} \right | - \gamma_{j} \right ]_
#' {+}^{2}}}}
#' and for the l0 penalty
#' \denq{\phi_{l_{0},m}(\gamma) = \max_{x\in S^{p}}\sum_{j=1}^{m}{\sum_{i=1}^{n}
#' { \left [ \left ( \mu_{j} a_{i}^{T}x_{j} \right )^{2} - \gamma_{j} \right ]_
#' {+}}}}
#' All of these functions are optimized using gradient decent. In this
#' implementation, a relative penalty is implemented with rho. Gamma is
#' calculated using rho and the maximal norm value.
#'
#'
#' @param A Input matrix of size (p x n) with p < n.
#' @param k Number of components, 0 < k < p.
#' @param rho Relative sparsity weight factor of the optimization. Either a
#' vector of floats of size k or float which will be repeated k times. 0 < rho
#' < 1.
#' @param penalty Penalty type to use in the optimization. Either "l0" or "l1".
#' The default is "l1" since it performed best in experiments.
#' @param center Centers the data. Either TRUE or FALSE. Default is TRUE.
#' @param block Optimization method. If FALSE, the components are calculated
#' individually. If TRUE, all components are calculated at the same time.
#' Default is FALSE.
#' @param mu Mean to be applied to each component in the block. Either a vector
#' of float of size k or a float which will be repeated k times. Only used if
#' block is TRUE. Default is FALSE.
#' @param iter_max Maximum iterations when adjusting components with
#' gradient descent. Default is 1000.
#' @param epsilon Epsilon of the gradient descent stopping function. Default is
#' 1e-4.
#'
#' @return List containing: \describe{
#'   \item{loadings}{The PCA components}
#'   \item{scores}{Scores of the components on A}
#'   \item{a_approx}{Reconstructed version of A using the components}
#'   \item{prop_sparse}{Proportion of sparsity of the components}
#'   \item{exp_var}{Explained ratio of variance of the components}
#'   \item{centers}{Centers of matrix A if center == TRUE}
#' }
#'
#' @examples
#' set.seed(360)
#' p <- 20
#' n <- 50
#' k <- 5
#' A <- scale(matrix(rnorm(p * n), nrow = p, ncol = n), scale = FALSE)
#' rho <- 0.1
#' # rho <- c(0.1, 0.2, 0.1, 0.2, 0.1)
#' mu <- 1
#' # mu <- c(1, 1.5, 0.5, 2, 1)
#'
#' # Single unit with l1 penalty
#' gpower(A, k, rho, 'l1', FALSE)
#'
#' # Single unit with l0 penalty
#' gpower(A, k, rho, 'l0', FALSE)
#'
#' # Block with l1 penalty
#' gpower(A, k, rho, 'l1', FALSE, TRUE, mu)
#'
#' # Block with l0 penalty
#' gpower(A, k, rho, 'l0', FALSE, TRUE, mu)
#'
#' @references
#' Journee, M., Nesterov, Y., Richtarik, P. and Sepulchre, R. (2010)
#' Generalized Power Method for Sparse Principal Component Analysis.
#' *Journal of Machine Learning Research. 11*, 517-553.
#'
#' @export
gpower <-
  function(A,
           k,
           rho,
           penalty = c("l0", "l1"),
           center = c(TRUE, FALSE),
           block = c(TRUE, FALSE),
           mu = NA,
           iter_max = 1000,
           epsilon = 1e-4) {
    UseMethod("gpower")
  }

#' @export
gpower.default <-
  function(A,
           k,
           rho,
           penalty = "l1",
           center = TRUE,
           block = FALSE,
           mu = NA,
           iter_max = 1000,
           epsilon = 1e-4) {
    A <- as.matrix(A)


    # Checks ------------------------------------------------------------------

    if (any(is.na(X))) {
      warning("Missing values are omitted: na.omit(X).")
      X <- stats::na.omit(X)
    }

    if (any(is.na(rho))) {
      warning("NA found in rho")
    }

    if (block == 1 & is.na(mu)) {
      warning("Mu needs to be defined for Block algoritm")
    }

    # Initialize gpower object ------------------------------------------------
    gpowerObj <- list(
      loadings = NULL,
      scores = NULL
    )

    p <- nrow(A)
    n <- ncol(A)

    iter_max <- iter_max
    epsilon <- epsilon

    Z <- matrix(0, nrow = n, ncol = k)

    if (length(rho) == 1) {
      rho <- rep(rho, k)
    }

    if (length(mu) == 1 & !is.na(mu)) {
      mu <- rep(mu, k)
    }
    # Add centering

    if (center) {
      A <- scale(A, scale=FALSE)
      gpowerObj$centers <- attr(A,"scaled:center")
    }

    X <- A

    # Single unit algorithm ---------------------------------------------------

    if (k == 1 | (k > 1 & !block)) {

      if (penalty == "l1") {

        for (c in 1:k) {
          # Loop on the components
          rho_c <- rho[c]
          norm_a_i <- rep(0, n)

          for (i in 1:n) {
            norm_a_i[i] <- norm(X[, i], "2")
          }

          rho_max <- max(norm_a_i)
          i_max <- which.max(norm_a_i)
          rho_c <- rho_c * rho_max

          # Initialisation
          x <- X[, i_max] / norm_a_i[i_max]
          f <- rep(0, iter_max)
          iter <- 1

          while (TRUE) {
            X_x <- t(X) %*% x
            t_resh <- sign(X_x) * pmax(abs(X_x) - rho_c, 0)

            # Cost function
            f[iter] <- sum(t_resh**2)

            if (f[iter] == 0) {
              # Sparsity is too high
              warning("Sparcity is set to high, all entries of loading vector are zero")
              break
            }
            else {
              gradient <- X %*% t_resh
              x <- gradient / norm(gradient, "2")
            }

            if (iter > 2) {
              # Stopping criteria
              if ((f[iter] - f[iter - 1]) / f[iter - 1] < epsilon |
                iter > iter_max) {
                if (iter > iter_max) {
                  print("Max iterations reached")
                  break
                }
                break
              }
            }

            iter <- iter + 1
          }

          X_x <- t(X) %*% x
          pattern <- (abs(X_x) - rho_c) > 0 # Pattern of sparsity
          z <- sign(X_x) * max(abs(X_x) - rho_c, 0)

          if (max(abs(z > 0))) {
            z <- z / norm(z, "2")
          }

          z <- pattern_filling(X, pattern, z)
          y <- X %*% z
          X <- X - y %*% z
          Z[, c] <- z
        }
      }

      if (penalty == "l0") {

        for (c in 1:k) {
          # Loop on the components
          rho_c <- rho[c]
          norm_a_i <- rep(0, n)

          for (i in 1:n) {
            norm_a_i[i] <- norm(X[, i], "2")
          }

          rho_max <- max(norm_a_i)
          i_max <- which.max(norm_a_i)
          rho_c <- rho_c * rho_max^2

          # Initialisation
          x <- X[, i_max] / norm_a_i[i_max]
          f <- rep(0, iter_max)
          iter <- 1

          while (TRUE) {
            X_x <- t(X) %*% x
            t_resh <- pmax(X_x^2 - rho_c, 0)

            # Cost function
            f[iter] <- sum(t_resh)

            if (f[iter] == 0) {
              # Sparsity is too high
              warning("Sparcity is set to high, all entries of loading vector are zero")
              break
            }
            else {
              gradient <- X %*% ((t_resh > 0) * X_x)
              x <- gradient / norm(gradient, "2")
            }

            if (iter > 2) {
              # Stopping criteria
              if ((f[iter] - f[iter - 1]) / f[iter - 1] < epsilon |
                iter > iter_max) {
                if (iter > iter_max) {
                  print("Max iterations reached")
                  break
                }
                break
              }
            }

            iter <- iter + 1
          }

          pattern <-
            ((t(X) %*% x)^2 - rho_c > 0) # Pattern of sparsity
          y <- x
          pattern_inv <- pattern == 0

          z <- t(X) %*% y
          z[pattern_inv] <- 0
          norm_z <- norm(z, "2")
          z <- z / norm_z
          y <- y * norm_z

          # Deflate
          X <- X - y %*% t(z)
          Z[, c] <- z
        }
      }
    }

    # Block algorithm with mu = 1 ---------------------------------------------
    else if (k > 1 & block & sum(mu == 1) == length(mu)) {
      norm_a_i <- rep(0, n)

      for (i in 1:n) {
        norm_a_i[i] <- norm(X[, i], "2")
      }

      i_max <- which.max(norm_a_i)

      # Initialization

      qr_decomp <- qr(cbind(
        X[, i_max] / norm_a_i[i_max],
        matrix(rnorm(p * (k - 1)), nrow = p)
      ),
      LAPACK = T
      ) # LAPACK to get MatLab qr(X, 0) results

      x <- qr.Q(qr_decomp)
      rho_max <- qr.R(qr_decomp)


      f <- rep(0, iter_max)
      iter <- 1

      if (penalty == "l1") {
        rho <- rho %*% rho_max

        while (TRUE) {
          X_x <- t(X) %*% x
          t_resh <- pmax(abs(X_x) - kronecker(matrix(1, n, 1), rho))
          f[iter] <- sum(t_resh^2)

          if (f[iter] == 0) {
            warning("Sparcity is set to high, all entries of loading vector are zero")
            break
          } else {
            gradient <- matrix(rep(0, p * k), nrow = p)

            for (i in 1:k) {
              pattern <- t(t_resh[, i]) > 0
              gradient[, i] <- X[, pattern] %*% (t_resh[pattern, i] * sign(X_x[pattern, i]))
            }

            svd_decomp <- svd(gradient)
            x <- svd_decomp$u %*% t(svd_decomp$v)
          }

          if (iter > 2) {
            if (iter > iter_max) {
              print("Max iterations reached")
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

        pattern <- abs(X_x) - kronecker(matrix(1, n, 1), rho) > 0
        Z <- pattern_filling(X, pattern, Z, mu)
      }

      if (penalty == "l0") {
        rho <- rho %*% (rho_max^2)

        while (TRUE) {
          X_x <- t(X) %*% x
          t_resh <- pmax(X_x^2 - kronecker(matrix(1, n, 1), rho))
          f[iter] <- sum(sum(t_resh))
          if (f[iter] == 0) {
            warning("Sparcity is set to high, all entries of loading vector are zero")
            break
          }
          else {
            gradient <- matrix(rep(0, p * k), nrow = p)

            for (i in 1:k) {
              pattern <- t(t_resh[, i]) > 0
              gradient[, i] <- X[, pattern] %*% X_x[pattern, i]
            }

            svd_decomp <- svd(gradient)
            x <- svd_decomp$u %*% t(svd_decomp$v)
          }

          if (iter > 2) {
            if (iter > iter_max) {
              print("Max iterations reached")
              break
            }
            if ((f[iter] - f[iter - 1]) / f[iter - 1] < epsilon) {
              break
            }
          }
          iter <- iter + 1
        }

        pattern <- (X_x^2 - kronecker(matrix(1, n, 1), rho)) > 0
        pattern_inv <- pattern == 0
        Z <- X_x
        Z[pattern_inv] <- 0
        norm_z <- rep(0, k)

        for (i in 1:k) {
          norm_z[i] <- norm(Z[, i], "2")
          if (norm_z[i] > 0) {
            Z[, i] <- Z[, i] / norm_z[i]
          }
        }
      }
    }

    # Block algorithm with mu != 1 ---------------------------------------------
    else if (k > 1 & block & sum(mu == 1) < length(mu)) {
      norm_a_i <- rep(0, n)

      for (i in 1:n) {
        norm_a_i[i] <- norm(X[, i], "2")
      }

      i_max <- which.max(norm_a_i)

      # Initialization

      qr_decomp <- qr(cbind(
        X[, i_max] / norm_a_i[i_max],
        matrix(rnorm(p * (k - 1)), nrow = p)
      ),
      LAPACK = T
      ) # LAPACK to get MatLab qr(X, 0) results

      x <- qr.Q(qr_decomp)
      rho_max <- qr.R(qr_decomp)


      f <- rep(0, iter_max)
      iter <- 1

      if (penalty == "l1") {
        rho <- rho * mu %*% rho_max

        while (TRUE) {
          X_x <- t(X) %*% x

          for (i in 1:k) {
            X_x[, i] <- X_x[, i] * mu[i]
          }

          t_resh <- pmax(abs(X_x) - kronecker(matrix(1, n, 1), rho))
          f[iter] <- sum(t_resh^2)

          if (f[iter] == 0) {
            warning("Sparcity is set to high, all entries of loading vector are zero")
            break
          } else {
            gradient <- matrix(rep(0, p * k), nrow = p)

            for (i in 1:k) {
              pattern <- t(t_resh[, i]) > 0
              gradient[, i] <- X[, pattern] %*% (t_resh[pattern, i] * sign(X_x[pattern, i]))
            }

            svd_decomp <- svd(gradient)
            x <- svd_decomp$u %*% t(svd_decomp$v)
          }

          if (iter > 2) {
            if (iter > iter_max) {
              print("Max iterations reached")
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
          X_x[,i] <- X_x[, i]*mu[i]
          Z[, i] <- sign(X_x[, i]) * max(abs(X_x[, i]) - rho[i], 0)
          if (max(abs(Z[, i]) > 0) > 0) {
            Z[, i] <- Z[, i] / norm(Z[, i], "2")
          }
        }

        pattern <- abs(X_x) - kronecker(matrix(1, n, 1), rho) > 0
        Z <- pattern_filling(X, pattern, Z, mu)
      }

      if (penalty == "l0") {
        rho <- rho * (mu %*% rho_max)^2

        while (TRUE) {
          X_x <- t(X) %*% x

          for (i in 1:k) {
            X_x[, i] <- X_x[, i] * mu[i]
          }

          t_resh <- pmax(X_x^2 - kronecker(matrix(1, n, 1), rho))
          f[iter] <- sum(sum(t_resh))
          if (f[iter] == 0) {
            warning("Sparcity is set to high, all entries of loading vector are zero")
            break
          }
          else {
            gradient <- matrix(rep(0, p * k), nrow = p)

            for (i in 1:k) {
              pattern <- t(t_resh[, i]) > 0
              gradient[, i] <- X[, pattern] %*% X_x[pattern, i]
            }

            svd_decomp <- svd(gradient)
            x <- svd_decomp$u %*% t(svd_decomp$v)
          }

          if (iter > 2) {
            if (iter > iter_max) {
              print("Max iterations reached")
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

        pattern <- (X_x^2 - kronecker(matrix(1, n, 1), rho)) > 0
        pattern_inv <- pattern == 0
        Z <- t(X) %*% x
        Z[pattern_inv] <- 0
        norm_z <- rep(0, k)

        for (i in 1:k) {
          norm_z[i] <- norm(Z[, i], "2")
          if (norm_z[i] > 0) {
            Z[, i] <- Z[, i] / norm_z[i]
          }
        }
      }
    }



    # Update gpower object ----------------------------------------------------
    weights <- Z # W
    scores <- A %*% weights # T
    P <- t(A) %*% ginv(t(scores))
    A_approx <- scores %*% t(P)

    gpowerObj$loadings <- Z
    gpowerObj$scores <- scores
    gpowerObj$a_approx <- A_approx
    gpowerObj$prop_sparse <- sum(rowSums(Z == 0)) / (n * k)

    # Explained variance ratio
    Svar <-
      1 - (norm((A - A_approx), type = "F") / norm(A, type = "F"))^2

    gpowerObj$exp_var <- Svar


    class(gpowerObj) <- "gpower"

    gpowerObj
  }


#' @export
print.gpower <- function(x, ...) {

  # Print gpower --------------------------------------------------------------

  cat("Proportion of Explained Variance\n")
  print(round(x$exp_var, 3))
  cat("\nSparse loadings:\n")
  print(round(x$loadings, 3))

}


#' @export
summary.gpower <- function(object, ...) {

  # Summary gpower ------------------------------------------------------------


  x <- t(data.frame(
    var = round(object$exp_var, 3),
    prop_sparse = round(object$prop_sparse, 3)
  ))

  rownames(x) <- c(
    "Proportion of Explained Variance",
    "Proportion of Sparsity"
  )

  x
}

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

    Z <-
      gpower(X, k, gamma, penalty, block, mu, iter_max, epsilon)

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

X <- as.matrix(read.csv("tests/testthat/data.csv", header = FALSE))


gpower(A=X, k = 5, rho = 0.1,  penalty = "l1")
gpower(A=X, k = 5, rho = 0.01, penalty = "l0")

gpower(A=X, k = 5, rho = 0.1, penalty = "l1", center=TRUE, block=TRUE, mu=1)
gpower(A=X, k = 5, rho = 0.01, penalty = "l0", center=TRUE, block=TRUE, mu=1)

gpower(A=X, k = 5, rho = 0.1, penalty = "l1", center=TRUE, block=TRUE,
                           mu=c(1,0.5,0.33,0.25,0.2))
gpower(A=X, k = 5, rho = 0.01, penalty = "l0", center=TRUE, block=TRUE,
                           mu=c(1,0.5,0.33,0.25,0.2))

# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#
# Code layout taken from the sparse-pca package by Benjamin Erichson


library(MASS)

#' @export
gpower <-
  function(A,
           k,
           rho,
           penalty = c('l0', 'l1'),
           block = c(TRUE, FALSE),
           mu = NA,
           iter_max = 1000,
           epsilon = 1e-4)
    UseMethod("gpower")

#' @export
gpower.default <-
  function(A,
           k = 1,
           rho = 0,
           penalty = 'l1',
           block = FALSE,
           mu = NA,
           iter_max = 1000,
           epsilon = 1e-4) {
    A <- as.matrix(A)
    X <- A

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Checks
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (any(is.na(X))) {
      warning("Missing values are omitted: na.omit(X).")
      X <- stats::na.omit(X)
    }

    if (block == 1 & is.na(mu)) {
      warning("Mu needs to be defined for Block algoritm")
      X <- stats::na.omit(X)
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Init gpower object
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gpowerObj = list(loadings = NULL,
                     scores = NULL,
                     testing = FALSE)

    p <- nrow(X)
    n <- ncol(X)

    iter_max <- iter_max
    epsilon <- epsilon

    Z <- matrix(0, nrow = n, ncol = k)

    if (length(rho) == 1) {
      rho <- rep(rho, k)
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Init Single unit algorithm
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (k == 1 | (k > 1 & !block)) {
      if (penalty == 'l1') {
        # L1 penalty

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
            f[iter] <- sum(t_resh ** 2)

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

        gpowerObj$testing <- TRUE
      }

      if (penalty == "l0") {
        warning("L0 penalty is not implemented yet")
      }
    }


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Update gpower object
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    weights <- Z # W
    scores <- A %*% weights # T
    P <- t(A) %*% ginv(t(scores))
    A_approx <- scores %*% t(P)

    Svar <-
      1 - (norm((A - A_approx), type = "F") / norm(A, type = "F")) ^ 2

    gpowerObj$loadings <- Z
    gpowerObj$scores <- scores
    gpowerObj$a_approx <- A_approx
    gpowerObj$prop_sparse <-  sum(rowSums(Z == 0)) / (n * k)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Explained variance and explained variance ratio
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    gpowerObj$exp_var <- Svar


    class(gpowerObj) <- "gpower"
    gpowerObj
  }


#' @export
print.gpower <- function(x , ...) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Print gpower
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  cat("Proportion of Explained Variance\n")
  print(round(x$exp_var, 3))
  cat("\nSparse loadings:\n")
  print(round(x$loadings, 3))
}


#' @export
summary.gpower <- function(object , ...)
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Summary gpower
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  x <- t(data.frame(
    var = round(object$exp_var, 3),
    prop_sparse = round(object$prop_sparse, 3)
  ))

  rownames(x) <- c('Explained variance',
                   'Proportion of Sparsity')

  # x <- as.matrix(x)

  x
}

#' @export
auto_gpower <- function(X,
                        k,
                        prop_sparse,
                        penalty = c('l0', 'l1'),
                        block = c(TRUE, FALSE),
                        mu = NA,
                        iter_max = 1000,
                        epsilon = 1e-4)
  UseMethod("auto_gpower")


#' @export
auto_gpower.default <- function(X,
                                k,
                                prop_sparse,
                                penalty = 'l1',
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
    upper <- max(norm(X, type = "2"))
  }
  if (penalty == "l1") {
    upper <- max(norm(X, type = "2")) ^ 2
  }
  i <- 0

  while (n_zeros_X != n_zeros_sparse & iter_max > i) {
    cut <- (lower + upper) / 2

    if (penalty == "l0") {
      gamma <- (cut ^ 2)
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

# Support function for gpower
pattern_filling <- function(A,
                            pattern,
                            Z = NA,
                            mu = NA) {
  # Compute a local maximizer of
  # max_{X,Z} trace(X^T A Z N)  s.t. X^T X=I_m and Z(P)=0 and Diag(Z^T Z)=I_m

  p <- nrow(A)
  n <- ncol(A)
  m <- ncol(pattern)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Single unit case
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (m == 1) {
    support <- pattern != 0

    if (sum(support) == 0) {
      z_red <- rep(0, n)
      support <- 1:n

    } else if (sum(support) == 1) {
      z_red <- 1

    } else {
      u <- Z[support]
      epsilon <- 1e-6
      iter_max <- 1000
      f <- rep(0, iter_max)
      iter <- 1
      A_red <- A[, support]

      while (TRUE) {
        temp <- t(A_red) %*% (A_red %*% u)
        u <- temp / norm(temp, "2")
        f[iter] <- -2 * t(u) %*% temp

        if (iter > 2) {
          # Stopping criteria
          if (abs(f[iter] - f[iter - 1]) / abs(f[iter - 1]) < epsilon |
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

      z_red <- u
    }

    z <- rep(0, n)
    z[support] <- z_red
    return(z)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # block case with mu_i == 1
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (is.na(mu) & !is.na(Z)) {
    # Invert pattern
    pattern_inv <- pattern - 1
    pattern_inv[pattern_inv == -1] <- 0
    iter_max <- 1000
    epsilon <- 1e-6
    f <- rep(0, iter_max)
    iter <- 1

    while (TRUE) {
      AZ <- A %*% Z
      # Creates d, u and v
      list2env(svd(AZ), .GlobalEnv)

      X <- u %*% t(v)
      ff <- 0
      for (i in 1:m) {
        ff <- ff + t(X[, i]) %*% AZ[, i]
      }

      f[iter] <- ff


      Z <- t(A) %*% X
      Z[pattern_inv] <- 0

      for (i in 1:m) {
        norm_Z <- norm(Z[, i], type = "2")
        if (norm_Z > 0) {
          Z[, i] <- Z[, 1] / norm_Z
        }
      }

      if (iter > 2) {
        if (iter > iter_max) {
          print("Maximum number of iterations reached")
          break
        }
        if (abs(f[iter] - f[iter - 1]) / abs(f[iter - 1]) < epsilon) {
          break
        }
      }

      iter <- iter + 1
    }
    return(Z)
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # General block case
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (!is.na(Z) & !is.na(mu)) {
    pattern_inv <- pattern - 1
    pattern_inv[pattern_inv == -1] <- 0
    iter_max <- 1000
    epsilon <- 1e-6
    f <- rep(0, iter_max)
    iter <- 1
    if (length(mu) == 1) {
      mu <- rep(mu, m)
    }

    while (TRUE) {
      AZ <- A %*% Z

      for (i in 1:m) {
        AZ[, i] <- AZ[, i] * mu[i]
      }

      # Creates d, u and v
      list2env(svd(AZ), .GlobalEnv)
      X <- u %*% t(v)

      ff <- 0
      for (i in 1:m) {
        ff <- ff + t(X[, i]) %*% AZ[, i]
      }
      f[iter] <- ff

      Z <- t(A) %*% X
      for (i in 1:m) {
        Z[, i] <- Z[, i] * mu[i]
      }
      Z[pattern_inv] <- 0

      for (i in 1:m) {
        norm_Z <- norm(Z[, 1], type = "2")
        if (norm_Z > 0) {
          Z[, i] <- Z[, i] / norm_Z
        }
      }

      if (iter > 2) {
        if (iter > iter_max) {
          print("Maximum number of iterations reached")
          break
        }
        if (abs(f[iter] - f[iter - 1]) / abs(f[iter - 1]) < epsilon) {
          break
        }
      }

      iter <- iter + 1
    }
    return(Z)
  }

}

library(sparsepca)

X <- as.matrix(read.csv("tests/data.csv", header = FALSE))


## Comparison
gpow <- gpower(X, rho = 0.1, k = 5)
sparse <- spca(X,
               k = 5,
               scale = FALSE,
               center = FALSE)

n_zeros_sparse <- sum(rowSums(gpow$loadings == 0))

## Auto
# auto <- auto_gpower(X, k=5, prop_sparse=0.8)
# summary(auto)
# print(auto)

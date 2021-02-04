# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#
# Code layout taken from the sparse-pca package by Benjamin Erichson


#' @export
gpower <-
  function(X,
           rho = 0,
           k = 1,
           penalty = 'l1',
           block = 0,
           mu = NA,
           iter_max = 1000,
           epsilon = 1e-5)
    UseMethod("gpower")

#' @export
gpower.default <-
  function(X,
           rho = 0,
           k = 1,
           penalty = 'l1',
           block = 0,
           mu = NA,
           iter_max = 1000,
           epsilon = 1e-4) {
    X <- as.matrix(X)

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
    gpowerObj = list(
      loadings = NULL,
      transform = NULL,
      scores = NULL,
      eigenvalues = NULL,
      testing = FALSE
    )

    p <- nrow(X)
    n <- ncol(X)

    iter_max <- iter_max
    epsilon <- epsilon

    Z <- matrix(0, nrow = n, ncol = k)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Init Single unit algorithm
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (k == 1 | (k > 1 & block == 0)) {
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

          # Initialisation
          x <- X[, i_max] / norm_a_i[i_max]
          f <- rep(0, iter_max)
          iter <- 1

          while (TRUE) {
            X_x <- t(X) %*% x
            t_resh <- sign(X_x) * max(abs(X_x) - rho_c, 0)

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

    gpowerObj$loadings <- Z
    # gpowerObj$transform <- A
    # gpowerObj$scores <- X %*% B
    # gpowerObj$eigenvalues <- svd_update$d / (n - 1)
    # gpowerObj$objective <- obj


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Explained variance and explained variance ratio
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # gpowerObj$sdev <-  sqrt( gpowerObj$eigenvalues )
    # gpowerObj$var <- sum( apply( Re(X) , 2, stats::var ) )
    # if(is.complex(X)) gpowerObj$var <- Re(gpowerObj$var + sum( apply( Im(X) , 2, stats::var ) ))


    class(gpowerObj) <- "gpower"
    return(gpowerObj)

  }


#' @export
print.gpower <- function(x , ...) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Print gpower
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # cat("Standard deviations:\n")
  # print(round(x$sdev, 3))
  # cat("\nEigenvalues:\n")
  # print(round(x$eigenvalues, 3))
  # cat("\nSparse loadings:\n")
  # print(round(x$loadings, 3))
  print(x$loadings)
}


#' @export
summary.gpower <- function(object , ...)
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Summary gpower
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # variance = object$sdev**2
  # explained_variance_ratio = variance / object$var
  # cum_explained_variance_ratio = cumsum( explained_variance_ratio )
  #
  # x <- t(data.frame( var = round(variance, 3),
  #                    sdev = round(object$sdev, 3),
  #                    prob = round(explained_variance_ratio, 3),
  #                    cum = round(cum_explained_variance_ratio, 3)))
  #
  # rownames( x ) <- c( 'Explained variance',
  #                     'Standard deviations',
  #                     'Proportion of variance',
  #                     'Cumulative proportion')
  #
  # colnames( x ) <- paste(rep('PC', length(object$sdev)), 1:length(object$sdev), sep = "")
  #
  # x <- as.matrix(x)

  return(object$loadings)
}

pattern_filling <- function(X,
                            pattern,
                            Z,
                            mu=NA) {

  # Compute a local maximizer of
  # max_{X,Z} trace(X^T A Z N)  s.t. X^T X=I_m and Z(P)=0 and Diag(Z^T Z)=I_m

  p <- nrow(X)
  n <- ncol(X)
  m <- ncol(pattern)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Single unit case
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (m == 1) {
    support <- pattern != 0

    if (sum(support) == 0) {
      z_red <- rep(0, n)
      support <- 1:n

    } else if (sum(support) == 1){
      z_red <- 1

    } else {
      u <- Z[support]
      epsilon <- 1e-6
      iter_max <- 1000
      f <- rep(0, iter_max)
      iter <- 1
      X_red <- X[, support]

      while (TRUE) {
        temp <- t(X_red) %*% (X_red %*% u)
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

  # TODO: Add block case with mu_i == 1

  # TODO: Add general block case
}

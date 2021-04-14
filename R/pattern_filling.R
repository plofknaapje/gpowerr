pattern_filling <- function(A, pattern, Z = NA, mu = NA) {
    # Compute a local maximizer of max_{X,Z} trace(X^T A Z N) s.t. X^T X=I_m and Z(P)=0 and Diag(Z^T Z)=I_m

    p <- nrow(A)
    n <- ncol(A)
    m <- ncol(pattern)

    # Single unit case ----------------------------------------------------------
    if (m == 1) {
        support <- pattern != 0

        if (sum(support) == 0) {
            z_red <- rep(0, n)
            support <- 1:n
        } else if (sum(support) == 1) {
            z_red <- 1
        } else {
            u <- Z[support]
            epsilon <- 1e-06
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

    # block case with mu_i == 1 -------------------------------------------------
    if (sum(is.na(mu)) == 0 & sum(mu == 1) == length(mu)) {
        # Invert pattern
        pattern_inv <- pattern == 0
        iter_max <- 1000
        epsilon <- 1e-06
        f <- rep(0, iter_max)
        iter <- 1

        while (TRUE) {
            AZ <- A %*% Z
            # Creates d, u and v
            svd_lst <- svd(AZ)
            u <- svd_lst$u
            v <- svd_lst$v

            X <- u %*% t(v)
            ff <- 0
            for (i in 1:m) {
                ff <- ff + t(X[, i]) %*% AZ[, i]
            }

            f[iter] <- ff

            Z <- t(A) %*% X
            Z[pattern_inv] <- 0

            for (i in 1:m) {
                norm_Z <- norm(Z[, i], "2")
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

    # General block case --------------------------------------------------------
    if (sum(is.na(mu)) == 0) {
        pattern_inv <- pattern - 1
        pattern_inv[pattern_inv == -1] <- 0
        iter_max <- 1000
        epsilon <- 1e-06
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
                norm_Z <- norm(Z[, 1], "2")
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

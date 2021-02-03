# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#
# Code layout taken from the sparse-pca package by Benjamin Erichson

#' @export
gpower <- function(X) UseMethod("gpower")
# spca <- function(X, k=NULL, alpha=1e-4, beta=1e-4, center=TRUE, scale=FALSE, max_iter=1000, tol=1e-5, verbose=TRUE)

#' @export
gpower.default <- function(X) {

  X <- as.matrix(X)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Checks
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (any(is.na(X))) {
    warning("Missing values are omitted: na.omit(X).")
    X <- stats::na.omit(X)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Init rpca object
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gpowerObj = list(loadings = NULL,
                 transform = NULL,
                 scores = NULL,
                 eigenvalues = NULL,
                 center = center,
                 scale = scale)

  n <- nrow(X)
  p <- ncol(X)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set target rank
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(k)) k <- min(n,p)
  if(k > min(n,p)) k <- min(n,p)
  if(k<1) stop("Target rank is not valid!")



  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Center/Scale data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(center == TRUE) {
    gpowerObj$center <- colMeans(X)
    X <- sweep(X, MARGIN = 2, STATS = gpowerObj$center, FUN = "-", check.margin = TRUE)
  } else { gpowerObj$center <- FALSE }

  if(scale == TRUE) {
    gpowerObj$scale <- sqrt(colSums(X**2) / (n-1))
    if(is.complex(gpowerObj$scale)) { gpowerObj$scale[Re(gpowerObj$scale) < 1e-8 ] <- 1+0i
    } else {gpowerObj$scale[gpowerObj$scale < 1e-8] <- 1}
    X <- sweep(X, MARGIN = 2, STATS = gpowerObj$scale, FUN = "/", check.margin = TRUE)
  } else { gpowerObj$scale <- FALSE }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute SVD for initialization of the Variable Projection Solver
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  svd_init <- svd(X)

  Dmax <- svd_init$d[1] # l2 norm

  A <- svd_init$v[,1:k]
  B <- svd_init$v[,1:k]

  V <- svd_init$v
  VD = sweep(V, MARGIN = 2, STATS = svd_init$d, FUN = "*", check.margin = TRUE)
  VD2 = sweep(V, MARGIN = 2, STATS = svd_init$d**2, FUN = "*", check.margin = TRUE)


  #--------------------------------------------------------------------
  #   Set Tuning Parameters
  #--------------------------------------------------------------------
  alpha <- alpha *  Dmax**2
  beta <- beta * Dmax**2

  nu <- 1.0 / (Dmax**2 + beta)
  kappa <- nu * alpha

  obj <- c()
  improvement <- Inf

  #--------------------------------------------------------------------
  #   Apply Variable Projection Solver
  #--------------------------------------------------------------------
  noi <- 1
  while (noi <= max_iter && improvement > tol) {

    # Update A:  X'XB = UDV'
    Z <- VD2 %*% (t(V) %*% B)
    svd_update <- svd(Z)
    A <- svd_update$u %*% t(svd_update$v)


    # Proximal Gradient Descent to Update B:
    grad <- VD2 %*% (t(V) %*% (A - B)) - beta * B
    B_temp <- B + nu * grad

    # l1 soft-threshold
    idxH <- which(B_temp > kappa)
    idxL <- which(B_temp <= -kappa)

    B = matrix(0, nrow = nrow(B_temp), ncol = ncol(B_temp))
    B[idxH] <- B_temp[idxH] - kappa
    B[idxL] <- B_temp[idxL] + kappa



    # compute residual
    R <- t(VD) - (t(VD) %*% B) %*% t(A)

    # compute objective function
    obj <- c(obj, 0.5 * sum(R**2) + alpha * sum(abs(B)) + 0.5 * beta * sum(B**2))

    # Break if obj is not improving anymore
    if(noi > 1){
      improvement <- (obj[noi-1] - obj[noi]) / obj[noi]
    }

    # Trace
    if(verbose > 0 && (noi-1) %% 10 == 0) {
      print(sprintf("Iteration: %4d, Objective: %1.5e, Relative improvement %1.5e", noi, obj[noi], improvement))
    }


    # Next iter
    noi <- noi + 1


  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Update gpower object
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gpowerObj$loadings <- B
  gpowerObj$transform <- A
  gpowerObj$scores <- X %*% B
  gpowerObj$eigenvalues <- svd_update$d / (n - 1)
  gpowerObj$objective <- obj


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Explained variance and explained variance ratio
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gpowerObj$sdev <-  sqrt( gpowerObj$eigenvalues )
  gpowerObj$var <- sum( apply( Re(X) , 2, stats::var ) )
  if(is.complex(X)) gpowerObj$var <- Re(gpowerObj$var + sum( apply( Im(X) , 2, stats::var ) ))


  class(gpowerObj) <- "gpower"
  return( gpowerObj )

}


#' @export
print.gpower <- function(x , ...) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Print rpca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("Standard deviations:\n")
  print(round(x$sdev, 3))
  cat("\nEigenvalues:\n")
  print(round(x$eigenvalues, 3))
  cat("\nSparse loadings:\n")
  print(round(x$loadings, 3))
}


#' @export
summary.gpower <- function( object , ... )
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Summary gpower
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  variance = object$sdev**2
  explained_variance_ratio = variance / object$var
  cum_explained_variance_ratio = cumsum( explained_variance_ratio )

  x <- t(data.frame( var = round(variance, 3),
                     sdev = round(object$sdev, 3),
                     prob = round(explained_variance_ratio, 3),
                     cum = round(cum_explained_variance_ratio, 3)))

  rownames( x ) <- c( 'Explained variance',
                      'Standard deviations',
                      'Proportion of variance',
                      'Cumulative proportion')

  colnames( x ) <- paste(rep('PC', length(object$sdev)), 1:length(object$sdev), sep = "")

  x <- as.matrix(x)

  return( x )
}

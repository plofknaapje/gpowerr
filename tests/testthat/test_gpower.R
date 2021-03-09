# devtools::use_package("testthat")

# Tests taken from the sparse-pca package by Benjamin Erichson

context("gpowerpca")

# Set seed
set.seed(1234)

# Generate Some Data in R -----------------------------------------------------
p <- 20
n <- 50
k <- 5
A <- scale(matrix(rnorm(p * n), nrow = p, ncol = n), scale = FALSE)
gamma <- rep(0.1, k)

# Test: G-Power PCA - penalty = l1, block = 0

# Test1: gpower returns matrix
testthat::test_that("Test 1: returns matrix", {
  testthat::expect_equal(is.matrix(gpower(A, rho = gamma, k = k, penalty = "l1", block = 0, mu = NA)$loadings), TRUE)
})

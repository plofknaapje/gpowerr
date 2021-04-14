# devtools::use_package('testthat')
context("gpowerr")

# Test1: gpower returns matrix
testthat::test_that("Test 1: Results hold up to matlab", {
    require(MASS)

    data <- as.matrix(read.csv("data.csv", header = FALSE))

    testthat::expect_gt(gpower(
        data = data,
        k = 5,
        rho = 0.1,
        penalty = "l1"
    )$exp_var,
    0.496)
    testthat::expect_gt(gpower(
        data = data,
        k = 5,
        rho = 0.01,
        penalty = "l0"
    )$exp_var,
    0.497)

    testthat::expect_gt(
        gpower(
            data = data,
            k = 5,
            rho = 0.1,
            penalty = "l1",
            center = TRUE,
            block = TRUE,
            mu = 1
        )$exp_var,
        0.49
    )
    testthat::expect_gt(
        gpower(
            data = data,
            k = 5,
            rho = 0.01,
            penalty = "l0",
            center = TRUE,
            block = TRUE,
            mu = 1
        )$exp_var,
        0.49
    )

    testthat::expect_gt(
        gpower(
            data = data,
            k = 5,
            rho = 0.1,
            penalty = "l1",
            center = TRUE,
            block = TRUE,
            mu = c(1, 0.5, 0.33, 0.25, 0.2)
        )$exp_var,
        0.45
    )
    testthat::expect_gt(
        gpower(
            data = data,
            k = 5,
            rho = 0.01,
            penalty = "l0",
            center = TRUE,
            block = TRUE,
            mu = c(1, 0.5, 0.33, 0.25, 0.2)
        )$exp_var,
        0.45
    )

    # Test for high values with blocks
    testthat::expect_true(
        is.list(gpower(
            data = data,
            k = 5,
            rho = 0.5,
            penalty = "l1",
            center = TRUE,
            block = TRUE,
            mu = 1
        ))
    )
    testthat::expect_true(
        is.list(gpower(
            data = data,
            k = 5,
            rho = 0.1,
            penalty = "l0",
            center = TRUE,
            block = TRUE,
            mu = 1
        ))
    )
})

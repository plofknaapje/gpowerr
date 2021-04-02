source("R/gpower.r")

gpower_var_plot <-
  function (data,
            k,
            intervals,
            penalty = "l1",
            center = TRUE,
            block = FALSE,
            mu = 1) {
    if (!block) {
      if (penalty == "l1") {
        values <- 0:(intervals - 1) / intervals
        objs <- lapply(values, function(x)
          gpower(
            data = data,
            k = k,
            rho = x,
            penalty = penalty,
            center = center
          ))
        exp_vars <- sapply(objs, function(x)
          x$exp_var)
        prop_sparses <- sapply(objs, function(x)
          x$prop_sparse)

      } else if (penalty == "l0") {
        values <- 0:(intervals - 1) / (intervals * 3)
        objs <- lapply(values, function(x)
          gpower(
            data = data,
            k = k,
            rho = x,
            penalty = penalty,
            center = center
          ))
        exp_vars <- sapply(objs, function(x)
          x$exp_var)
        prop_sparses <- sapply(objs, function(x)
          x$prop_sparse)
      } else {
        warning("penalty not recognized")
      }

    } else {
      if (penalty == "l1") {
        objs <- list()
        values <- 0:(intervals - 1) / intervals

        for (x in values) {
          print(x)
          obj <- tryCatch(
            gpower(
              data = data,
              k = k,
              rho = x,
              penalty = penalty,
              center = center,
              block = block,
              mu = mu
            ),
            error = function(e) {
              print(e)
              NULL
            }
          )
          if (is.list(obj)) {
            objs <- c(objs, list(obj))
          } else {
            print(obj)
            break
          }
        }
        exp_vars <- sapply(objs, function(x)
          x$exp_var)
        prop_sparses <- sapply(objs, function(x)
          x$prop_sparse)
        values <- values[1:length(exp_vars)]

      } else if (penalty == "l0") {
        objs <- list()
        values <- 0:(intervals - 1) / (intervals * 3)

        for (x in values) {
          print(x)
          obj <- tryCatch(
            gpower(
              data = data,
              k = k,
              rho = x,
              penalty = penalty,
              center = center,
              block = block,
              mu = mu
            ),
            error = function(e) {
              print(e)
              NULL
            }
          )
          if (is.list(obj)) {
            objs <- c(objs, list(obj))
          } else {
            print(obj)
            break
          }
        }
        exp_vars <- sapply(objs, function(x)
          x$exp_var)
        prop_sparses <- sapply(objs, function(x)
          x$prop_sparse)
        values <- values[1:length(exp_vars)]

      } else {
        warning("penalty not recognized")
      }
    }
    lbound <- floor(intervals/3)
    rbound <- 2*floor(intervals/3)

    alignment <- "topleft"
    if (sum(exp_vars[1:lbound] >= 0.8) +
        sum(prop_sparses[1:lbound] >= 0.8) >= 1) {
      alignment <- "bottomright"
      if (sum(exp_vars[rbound:length(exp_vars)] <= 0.2) +
          sum(prop_sparses[rbound:length(prop_sparses)] <= 0.2) >= 1) {
      alignment <- "right"
      }
    }

    plot(
      values,
      exp_vars,
      type = "l",
      main = "Variance and Sparseness as a function of Rho",
      xlab = "Rho",
      ylab = "Proportion",
      col = "blue",
      ylim = c(min(exp_vars, prop_sparses), 1)
    )
    lines(values, prop_sparses, type = "l", col = "orange")
    legend(
      x = alignment,
      legend = c("Explained Variance", "Sparseness"),
      col = c("blue", "orange"),
      lty = c(1, 1)
    )
  }

# gpower_var_plot(
#   data = read.csv("tests/testthat/data.csv", header = FALSE),
#   k = 15,
#   intervals = 40,
#   penalty = "l1",
#   center = TRUE,
#   block = FALSE,
#   mu = 1
# )

gpower_comp_var_plot <-
  function (data,
            k,
            intervals,
            penalty = "l1",
            center = TRUE,
            block = FALSE,
            mu = 1) {

    if (!block) {
      if (penalty == "l1") {
        values <- 0:(intervals - 1) / intervals
        objs <- lapply(values, function(x)
          gpower(
            data = data,
            k = k,
            rho = x,
            penalty = penalty,
            center = center
          ))
        comp_vars <-  t(sapply(objs, function(x) x$comp_var))
        print(comp_vars)

      } else if (penalty == "l0") {
        values <- 0:(intervals - 1) / (intervals * 3)
        objs <- lapply(values, function(x)
          gpower(
            data = data,
            k = k,
            rho = x,
            penalty = penalty,
            center = center
          ))
        comp_vars <- t(sapply(objs, function(x)
          x$comp_var))

      } else {
        warning("penalty not recognized")
      }

    } else {
      if (penalty == "l1") {
        objs <- list()
        values <- 0:(intervals - 1) / intervals

        for (x in values) {
          print(x)
          obj <- tryCatch(
            gpower(
              data = data,
              k = k,
              rho = x,
              penalty = penalty,
              center = center,
              block = block,
              mu = mu
            ),
            error = function(e) {
              print(e)
              NULL
            }
          )
          if (is.list(obj)) {
            objs <- c(objs, list(obj))
          } else {
            print(obj)
            break
          }
        }
        comp_vars <- t(sapply(objs, function(x)x$comp_var))
        values <- values[1:length(exp_vars)]

      } else if (penalty == "l0") {
        objs <- list()
        values <- 0:(intervals - 1) / (intervals * 3)

        for (x in values) {
          print(x)
          obj <- tryCatch(
            gpower(
              data = data,
              k = k,
              rho = x,
              penalty = penalty,
              center = center,
              block = block,
              mu = mu
            ),
            error = function(e) {
              print(e)
              NULL
            }
          )
          if (is.list(obj)) {
            objs <- c(objs, list(obj))
          } else {
            print(obj)
            break
          }
        }
        comp_vars <- t(sapply(objs, function(x)
          x$exp_var))

        values <- values[1:length(exp_vars)]

      } else {
        warning("penalty not recognized")
      }
    }
    print(values)
    print(comp_vars)

    matplot(
      values,
      comp_vars,
      type = "l",
      main = "Variance explained by each component",
      xlab = "Rho",
      ylab = "Variance"
    )
    legend('topright', legend=1:5,
           pch=1, horiz=TRUE, col=1:5)

  }


gpower_component_heatmap <- function(data,
                                     k,
                                     rho,
                                     penalty = "l1",
                                     center = TRUE,
                                     block = FALSE,
                                     mu = 1){
  pow <- gpower(data, k, rho, penalty, center, block, mu)
  heatmap(pow$loadings, Colv=NA)
}

# gpower_comp_var_plot(
#   data = read.csv("tests/testthat/data.csv", header = FALSE),
#   k = 5,
#   intervals = 40,
#   penalty = "l0",
#   center = TRUE,
#   block = FALSE,
#   mu = 1
# )


# gpower_component_heatmap(data = read.csv("tests/testthat/data.csv", header = FALSE),
#                            k = 5,
#                            rho = 0.1,
#                            penalty = "l1",
#                            center = TRUE,
#                            block = FALSE,
#                            mu = 1)

source("R/gpower.r")
library(qgraph)

data <- read.csv("tests/testthat/data.csv", header=FALSE)
penalty <- "l1"
block <- TRUE
mu <- 1
intervals <- 40
k <- 5
center <- TRUE

if (!block) {
  if (penalty == "l1") {
    values <- 0:(intervals-1)/intervals
    objs <- lapply(values, function(x)gpower(A=data, k=k, rho=x,
                                            penalty=penalty, center=center))
    exp_vars <- sapply(objs, function(x)x$exp_var)
    prop_sparses <- sapply(objs, function(x)x$prop_sparse)

  } else if (penalty == "l0") {
    values <- 0:(intervals-1)/(intervals*3)
    objs <- lapply(values, function(x)gpower(A=data, k=k, rho=x,
                                             penalty=penalty, center=center))
    exp_vars <- sapply(objs, function(x)x$exp_var)
    prop_sparses <- sapply(objs, function(x)x$prop_sparse)
  } else {
    warning("penalty not recognized")
  }

} else {
    if (penalty == "l1") {
      objs <- list()
      values <- 0:(intervals-1)/intervals

      for (x in values) {
        print(x )
        obj <- tryCatch(
          gpower(A=data, k=k, rho=x,
                   penalty=penalty, center=center, block=block, mu=mu),
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
      exp_vars <- sapply(objs, function(x)x$exp_var)
      prop_sparses <- sapply(objs, function(x)x$prop_sparse)
      values <- values[1:length(exp_vars)]

    } else if (penalty == "l0") {
      objs <- list()
      values <- 0:(intervals-1)/(intervals*3)

      for (x in values) {
        print(x )
        obj <- tryCatch(
          gpower(A=data, k=k, rho=x,
                 penalty=penalty, center=center, block=block, mu=mu),
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
      exp_vars <- sapply(objs, function(x)x$exp_var)
      prop_sparses <- sapply(objs, function(x)x$prop_sparse)
      values <- values[1:length(exp_vars)]

    } else {
      warning("penalty not recognized")
    }
}

plot(values, exp_vars, type="l", main="Variance and Sparseness as a function of Rho",
     xlab="Rho", ylab="Proportion", col="blue", ylim = c(min(exp_vars, prop_sparses),1))
lines(values, prop_sparses, type="l", col="orange")
legend(x="right", legend=c("Explained Variance", "Sparseness"), col=c("blue","orange"), lty=c(1,1))


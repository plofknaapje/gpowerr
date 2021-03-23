source("R/gpower.r")
library(qgraph)

data <- read.csv("tests/testthat/data.csv", header=FALSE)
values <- (1:19/20)

var <- sapply(values, function(x)gpower(A=data, k=6, rho=x,
                                        penalty="l1", center=TRUE)$exp_var)

plot(values, var, type="l", main="Explained Variance as a function of Rho",
     xlab="Rho", ylab="Explained Variance")
print(var)

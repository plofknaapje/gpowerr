# Matlab model test
library(gpowerpca)
library(sparsepca)

A <- read.csv("tests/data.csv", header = FALSE)

Z1 <- gpower(A, rho=0.1, k=5)$loadings
Z2 <- spca(A, k=5)

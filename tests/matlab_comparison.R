
X <- as.matrix(read.csv("tests/data.csv", header = FALSE))


## Comparison
gpow_single <- gpower(X,
  rho = 0.1,
  k = 5,
  penalty = "l1"
)

print(gpow_single)

## Auto
# auto <- auto_gpower(X, k = 5, prop_sparse = 0.4)
# summary(auto)

gpow <- gpower(X,
  rho = 0.1,
  k = 5,
  penalty = "l1",
  block = TRUE,
  mu = 1
)

print(gpow)

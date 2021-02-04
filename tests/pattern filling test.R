library(gpowerpca)
library(sparsepca)



m = 1000
V1 = rnorm(m, 0, 290)
V2 = rnorm(m, 0, 300)
V3 = -0.1*V1 + 0.1*V2 + rnorm(m,0,100)

X = cbind(V1,V1,V1,V1, V2,V2,V2,V2, V3,V3)
X = X + matrix(rnorm(length(X),0,1), ncol = ncol(X), nrow = nrow(X))
X = scale(X, scale = FALSE)

k <- 3
gamma <- rep(0.2, k)

gp_out<- gpower(X, rho=gamma, k=k)
spca_out <- spca(X, k=k, alpha=0, beta=0, center = TRUE, scale = FALSE, verbose=0)

gp_out$loadings
spca_out$loadings

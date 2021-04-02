source("R/var_plot.r")
library(qgraph)
data(big5)

colnames(big5)

pow <- gpower(data=big5, k=10, rho=0.01, penalty="l0", center=TRUE)

typeof(pow$loadings)

#gpower_var_plot(data=big5, k=5, intervals = 20, penalty = "l0")

gpower_component_heatmap(data=big5, k=5, rho=0.1, penalty="l1", center=TRUE)

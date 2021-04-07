source("R/plotting.r")
library(qgraph)
data(big5)

pow <- gpower(data = big5, k = 2, rho = 0.15, penalty = "l1", center = TRUE)

# gpower_var_plot(data = big5, k = 2, intervals = 20, penalty = 'l1')

gpower_comp_var_plot(data = big5, k = 5, intervals = 20, penalty = 'l1')

# gpower_component_heatmap(data = big5, k = 1, rho = 0.1, penalty = 'l1', center = TRUE)

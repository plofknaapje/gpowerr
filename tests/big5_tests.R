source("R/plotting.r")
library(qgraph)
data(big5)

# pow <- gpower(data = big5, k = 2, rho = 0.1, penalty = "l0", center = TRUE,
              block = TRUE, mu = 1)

# auto_gpower(data = big5, k = 5, prop_sparse = 0.2, penalty = "l0", block = TRUE, mu = 1)

# gpower_var_plot(data = big5, k =5, intervals = 20, penalty = 'l0')

gpower_comp_var_plot(data = big5, k = 5, intervals = 20, penalty = 'l1')

# gpower_component_heatmap(data = big5, k = 1, rho = 0.1, penalty = 'l1', center = TRUE)

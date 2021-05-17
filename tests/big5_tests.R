source("R/plotting.r")
library(qgraph)
data(big5)


pow <- gpower(data = big5, k = 2, rho = 0.1, penalty = "l0", center = TRUE,
              block = TRUE, mu = 1)

auto_gpower(data = big5, k = 5, prop_sparse = 0.5, penalty = "l1", block = FALSE, mu = 1)

gpower_var_plot(data = big5, k =5, intervals = 20, penalty = 'l1')

gpower_comp_var_plot(data = big5, k = 5, intervals = 20, penalty = 'l1')

row_labels = data.frame(label=substring(colnames(big5)[1:40],1,1), row.names = colnames(big5)[1:40])
gpower_component_heatmap(data = big5[,1:40], k = 5, rho = 0.01, penalty = 'l1', center = TRUE,
                         cluster_variables = TRUE, show_variable_names = TRUE,
                         variable_highlight = row_labels)

source("R/gpower.r")

#' Plots explained variance and sparsity of the components
#'
#' This plot shows what happens to the proportion of explained variance and
#' proportion of sparsity if rho increases. The sparsity is non-decreasing, but
#' the explained variance is not.
#' This is also the case in the original MatLab code of the inventers of the
#' method.
#'
#' @inheritParams gpower
#' @param intervals The amount of intervals in the range of values of rho for
#'   which gpower is run. For l1, the range of values is 0-1, for l0, the range
#'   is 0-0.33. The block algorithms tend to stop at 0.4 and will only cover the
#'   range where gpower is able to run.
#'
#' @examples
#' set.seed(360)
#' p <- 20
#' n <- 50
#' k <- 5
#' data <- scale(matrix(stats::rnorm(p * n), nrow = p, ncol = n), scale = FALSE)
#'
#' gpower_var_plot(
#'   data = data,
#'   k = k,
#'   intervals = 40,
#'   penalty = 'l1',
#'   center = TRUE,
#' )
#' @export
gpower_var_plot <-
    function(data,
             k,
             intervals,
             penalty = "l1",
             center = TRUE,
             block = FALSE,
             mu = 1) {
        if (!block) {
            if (penalty == "l1") {
                values <- 0:(intervals - 1) / intervals
                objs <-
                    lapply(values, function(x)
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
                objs <-
                    lapply(values, function(x)
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
                    obj <-
                        tryCatch(
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
                    obj <-
                        tryCatch(
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
        lbound <- floor(intervals / 3)
        rbound <- 2 * floor(intervals / 3)

        alignment <- "topleft"
        if (sum(exp_vars[1:lbound] >= 0.8) + sum(prop_sparses[1:lbound] >= 0.8) >= 1) {
            alignment <- "bottomright"
            if (sum(exp_vars[rbound:length(exp_vars)] <= 0.2) + sum(prop_sparses[rbound:length(prop_sparses)] <= 0.2) >= 1) {
                alignment <- "right"
            }
        }

        plot(
            values,
            exp_vars,
            type = "l",
            main = "Variance and Sparseness as a function of Rho",
            xlab = "Proportion of the upper bound",
            ylab = "Proportion",
            col = "blue",
            ylim = c(min(exp_vars, prop_sparses),
                     1)
        )
        graphics::lines(values, prop_sparses, type = "l", col = "orange")
        graphics::legend(
            x = alignment,
            legend = c("Explained Variance", "Sparseness"),
            col = c("blue", "orange"),
            lty = c(1, 1)
        )
    }

#' Plots explained variance of each component
#'
#' This plot shows what happens to the variance explained by each component if
#' rho is increased. The individual explained variances are neither
#' non-decreasing nor non-increasing. Their sum does trend downwards, but it is
#' also not nonincreasing just like the proportion of explained variance.
#'
#' @inheritParams gpower
#' @param intervals The amount of intervals in the range of values of rho for
#'   which gpower is run. For l1, the range of values is 0-1, for l0, the range
#'   is 0-0.33. The block algorithms tend to stop at 0.4 and will only cover the
#'   range where gpower is able to run.
#'
#' @examples
#' set.seed(360)
#' p <- 20
#' n <- 50
#' k <- 5
#' data <- scale(matrix(stats::rnorm(p * n), nrow = p, ncol = n), scale = FALSE)
#'
#' gpower_comp_var_plot(
#'   data = data,
#'   k = k,
#'   intervals = 40,
#'   penalty = 'l1',
#'   center = TRUE,
#' )
#' @export
gpower_comp_var_plot <-
    function(data,
             k,
             intervals,
             penalty = "l1",
             center = TRUE,
             block = FALSE,
             mu = 1) {
        if (!block) {
            if (penalty == "l1") {
                values <- 0:(intervals - 1) / intervals
                objs <-
                    lapply(values, function(x)
                        gpower(
                            data = data,
                            k = k,
                            rho = x,
                            penalty = penalty,
                            center = center
                        ))
                comp_vars <- t(sapply(objs, function(x)
                    x$comp_var))

            } else if (penalty == "l0") {
                values <- 0:(intervals - 1) / (intervals * 3)
                objs <-
                    lapply(values, function(x)
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
                    obj <-
                        tryCatch(
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
                    x$comp_var))
                values <- values[1:length(comp_vars)]

            } else if (penalty == "l0") {
                objs <- list()
                values <- 0:(intervals - 1) / (intervals * 3)

                for (x in values) {
                    obj <-
                        tryCatch(
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
                        break
                    }
                }
                comp_vars <- t(sapply(objs, function(x)
                    x$exp_var))

                values <- values[1:length(comp_vars)]

            } else {
                warning("penalty not recognized")
            }
        }

        graphics::matplot(
            values,
            comp_vars,
            type = "l",
            main = "Variance explained by each component",
            xlab = "Proportion of the upper bound",
            ylab = "Variance"
        )
        graphics::legend(
            "topright",
            legend = 1:k,
            pch = 1,
            horiz = TRUE,
            col = 1:5
        )

    }

#' Plots value of each column of the data for each component
#'
#' This plot shows the values of the components, connected to the columns of the
#' original data. The column data is ordered using builtin dendrogram of the
#' base heatmap function.
#'
#' @inheritParams gpower
#'
#' @examples
#' set.seed(360)
#' p <- 20
#' n <- 50
#' k <- 5
#' data <- scale(matrix(stats::rnorm(p * n), nrow = p, ncol = n), scale = FALSE)
#'
#' gpower_component_heatmap(
#'   data = data,
#'   k = 5,
#'   rho = 0.1,
#'   penalty = 'l1',
#'   center = TRUE,
#'   block = FALSE,
#'   mu = 1,
#'   cluster_variables = FALSE,
#'   show_variable_names = TRUE
#' )
#' @export
gpower_component_heatmap <-
    function(data,
             k,
             rho,
             penalty = "l1",
             center = TRUE,
             block = FALSE,
             mu = 1,
             variable_highlight = NA,
             cluster_variables = FALSE,
             show_variable_names = TRUE) {
        # Add custom colorbrew
        pow <- gpower(data, k, rho, penalty, center, block, mu)
        if (any(is.na(variable_highlight))){
            pheatmap::pheatmap(pow$loadings, cluster_cols = FALSE,
                               cluster_rows = cluster_variables,
                               show_rownames = show_variable_names,
                               main = "Heatmap of GPower component weights")
        } else {
            pheatmap::pheatmap(pow$loadings, cluster_cols = FALSE,
                               cluster_rows = cluster_variables,
                               show_rownames = show_variable_names,
                               annotation_row = variable_highlight,
                               main = "Heatmap of GPower component weights",
                               annotation_names_row = FALSE)
        }

        # Add legend!!!
    }

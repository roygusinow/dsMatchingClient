#' Plot Empirical Cumulative Distribution Functions (eCDFs) for Matched and Unmatched Data
#'
#' This function plots eCDFs of variables defined in a formula for both unmatched and matched datasets
#' using `ds.eCDF` from the DataSHIELD environment. It visually compares the distribution of covariates
#' between treated and control groups before and after matching.
#'
#' @param unmatched A character string. Name of the unmatched data frame on the server.
#' @param matched A character string. Name of the matched data frame on the server.
#' @param formula An R formula object specifying the variables. The first variable should be the treatment indicator.
#' @param len Integer. Number of interval points used to construct the domain of the eCDF.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function attempts to find active connections via `datashield.connections_find()`.
#'
#' @return A combined `ggplot2` object showing eCDF plots.
#'
#' @details
#' This function creates server-side subsets for treated and control groups using a binary treatment variable.
#' It then calls `ds.eCDF` to retrieve the cumulative distribution information and plots the results for each
#' covariate in the formula using `ggplot2`.
#'
#' @author Roy Gusinow
#'
#' @export
ds.eCDF_plot <- function(unmatched,
                      matched,
                      formula,
                      len = 40,
                      datasources = NULL) {

  check_formula(formula)

  f_vars <- all.vars(formula)
  treatment <- f_vars[1]

  plot_list <- c()
  for (df in c(unmatched, matched)){

    # get correct vectors for comparisons server side
    treatment_vec <- paste0(df, "$", treatment)
    assign_ZEROS(treatment_vec, datasources)
    assign_ONES(treatment_vec, datasources)

    # subset frames
    dsBaseClient::ds.dataFrameSubset(df.name = df,
                       V1.name = treatment_vec,
                       V2.name = 'ZEROS',
                       Boolean.operator = "==",
                       newobj = "control",
                       datasources = datasources)
    dsBaseClient::ds.dataFrameSubset(df.name = df,
                       V1.name = treatment_vec,
                       V2.name = 'ONES',
                       Boolean.operator = "==",
                       newobj = "treated",
                       datasources = datasources)

    # plot each variable in control and matched - changes because new objects on the server are created per loop
    one_plot_list <- lapply(f_vars[2:length(f_vars)], get_eCDF_plot_per_var, len = len, datasources = datasources)
    one_plot_list[[1]] <- one_plot_list[[1]] + ggplot2::ggtitle(ifelse(df == unmatched, "All", "Matched"))
    plot_list <- c(plot_list, one_plot_list)

  }

  # plot specifications
  grid_plot <- cowplot::plot_grid(plotlist = plot_list, byrow = FALSE, nrow = length(f_vars)-1, ncol = 2)

  x.grob <- grid::textGrob("Control Units",
                     gp=grid::gpar(fontface="bold", col="black", fontsize=15))
  y.grob <- grid::textGrob("Treated Units",
                     gp=grid::gpar(fontface="bold", col="black", fontsize=15), rot=90)
  title.grob <- grid::textGrob("eCDF Plot",
                         gp=grid::gpar(fontface="bold", col="black", fontsize=15))

  grid_plot <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grid_plot, left = y.grob, bottom = x.grob, top = title.grob))


  return(grid_plot)
}

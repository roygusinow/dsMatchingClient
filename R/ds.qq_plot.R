#' This function generates quantile-quantile (QQ) plots comparing the distributions of covariates
#' between treated and control units, both before and after matching.
#'
#' @param unmatched A character string. The name of the unmatched (pre-matching) data frame stored on the server.
#' @param matched A character string. The name of the matched data frame stored on the server.
#' @param formula A formula object. Specifies the treatment variable on the left-hand side, and covariates to compare on the right-hand side.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function uses `datashield.connections_find()` to detect active connections.
#'
#' @return A combined `ggplot` object consisting of side-by-side QQ plots for each covariate, comparing treated vs control units in the unmatched and matched datasets.
#'
#' @details
#' For each covariate, the function creates QQ plots of treated vs control distributions,
#' before and after matching. It subsets the data by treatment status server-side using
#' `ds.dataFrameSubset`, and compares covariate distributions using empirical quantiles
#' aggregated with `ds.quantile`.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item None, but uses `ds.dataFrameSubset` to create subsets of treated and control units.
#' }
#'
#' @author Roy Gusinow
#' @export
ds.qq_plot <- function(unmatched = "DST",
                       matched = "matched_pooled",
                       formula = NULL,
                       datasources = NULL){

  check_formula(formula)

  f_vars <- all.vars(formula)
  treatment <- f_vars[1]

  plot_list <- c()
  for (df in c(unmatched, matched)){

    # get correct vectors for comparisons server side
    treatment_vec <- paste0(df, "$", treatment)
    assign_ZEROS(treatment_vec, datasources = datasources)
    assign_ONES(treatment_vec, datasources = datasources)

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
    one_plot_list <- lapply(f_vars[2:length(f_vars)], get_qq_plot_per_var, datasources = datasources)
    one_plot_list[[1]] <- one_plot_list[[1]] + ggplot2::ggtitle(ifelse(df == unmatched, "All", "Matched"))
    plot_list <- c(plot_list, one_plot_list)

  }

  # plot specifications
  grid_plot <- cowplot::plot_grid(plotlist = plot_list, byrow = FALSE, nrow = length(f_vars)-1, ncol = 2)

  x.grob <- grid::textGrob("Control Units",
                     gp=grid::gpar(fontface="bold", col="black", fontsize=15))
  y.grob <- grid::textGrob("Treated Units",
                     gp=grid::gpar(fontface="bold", col="black", fontsize=15), rot=90)
  title.grob <- grid::textGrob("eQQ Plot",
                         gp=grid::gpar(fontface="bold", col="black", fontsize=15))

  grid_plot <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grid_plot, left = y.grob, bottom = x.grob, top = title.grob))


  return(grid_plot)
}


#' @name utils
#' @keywords internal

#' @aliases check_formula
#' @noRd
check_formula <- function(formula){
  # check if the formula is a valid formula object
  if (!inherits(formula, "formula")) {
    stop("The 'formula' argument must be a valid formula object.", call. = FALSE)
  }
  return(TRUE)
}

# useful functions
#' @noRd
assign_ZEROS <- function(vec = NULL, datasources){
  # assign zeros of the correct vec length to a server
  dsBaseClient::ds.make(toAssign = paste0(vec, "-", vec),
          newobj = "ZEROS", datasources = datasources)
}
#' @noRd
assign_ONES <- function(vec = NULL, datasources){
  # assign zeros of the correct vec length to a server
  dsBaseClient::ds.make(toAssign = paste0(vec, "-", vec, "+ 1"),
          newobj = "ONES", datasources = datasources)
}
#' @noRd
avg_compute_estimate <- function(fit,
                                 data,
                                 treatment = NULL,
                                 estimand = "ATE",
                                 comparison = "difference",
                                 datasources = NULL){
  # function to retrieve federated estimate

  dsBaseClient::ds.assign(toAssign = paste0(data, "$", treatment, " - ", data, "$", treatment),
            newobj = "ZEROS",
            datasources = datasources)
  if (estimand == "ATE"){
    data <- data
  }else if (estimand == "ATT"){

    dsBaseClient::ds.dataFrameSubset(df.name = data,
                       V1.name = paste0(data, "$", treatment),
                       V2.name = "ZEROS",
                       Boolean.operator = "!=",
                       keep.NAs = FALSE,
                       newobj = "new_data",
                       datasources = datasources, #all servers are used
                       notify.of.progress = FALSE) # remove treatment
    data <- "new_data"

  }else if (estimand == "ATC"){

    dsBaseClient::ds.dataFrameSubset(df.name = data,
                       V1.name = paste0(data, "$", treatment),
                       V2.name = "ZEROS",
                       Boolean.operator = "==",
                       keep.NAs = FALSE,
                       newobj = "new_data",
                       datasources = datasources, #all servers are used
                       notify.of.progress = FALSE) # remove treatment
    data <- "new_data"
  }

  # replace treatment with 0
  dsBaseClient::ds.assign(toAssign = data,
            newobj = "df_off",
            datasources = datasources)
  dsBaseClient::ds.assign(toAssign = paste0(data, "$", treatment, " - ", data, "$", treatment),
            newobj = treatment,
            datasources = datasources)
  dsBaseClient::ds.dataFrameSubset(df.name = "df_off",
                     V1.name = paste0(data, "$", treatment),
                     V2.name = paste0(data, "$", treatment),
                     Boolean.operator = "==",
                     rm.cols = which(dsBaseClient::ds.colnames("df_off", datasources = datasources)[[1]] == treatment),
                     keep.NAs = FALSE,
                     newobj = "df_off",
                     datasources = datasources, #all servers are used
                     notify.of.progress = FALSE) # remove treatment
  dsBaseClient::ds.cbind(c("df_off", treatment),
           newobj = "df_off",
           datasources = datasources)

  # replace treatment with 1
  dsBaseClient::ds.assign(toAssign = data,
            newobj = "df_on",
            datasources = datasources)
  dsBaseClient::ds.assign(toAssign = paste0(data, "$", treatment, " + !", data, "$", treatment),
            newobj = treatment,
            datasources = datasources)
  dsBaseClient::ds.dataFrameSubset(df.name = "df_on",
                     V1.name = paste0(data, "$", treatment),
                     V2.name = paste0(data, "$", treatment),
                     Boolean.operator = "==",
                     rm.cols = which(dsBaseClient::ds.colnames("df_on", datasources = datasources)[[1]] == treatment),
                     keep.NAs = FALSE,
                     newobj = "df_on",
                     datasources = datasources, #all servers are used
                     notify.of.progress = FALSE) # remove treatment
  dsBaseClient::ds.cbind(c("df_on", treatment),
           newobj = "df_on",
           datasources = datasources)


  coef_vec <- fit$coefficients[,1]
  names(coef_vec)[names(coef_vec) == "(Intercept)"] <- "Intercept"

  ds.genProp(formula = stats::as.formula(fit$formula), coefficients = coef_vec, data = "df_off", newobj = "df_off_pred", link = fit$family$link, datasources = datasources)
  ds.genProp(formula = stats::as.formula(fit$formula), coefficients = coef_vec, data = "df_on", newobj = "df_on_pred", link = fit$family$link, datasources = datasources)

  avg_off <- dsBaseClient::ds.mean("df_off_pred", type = "combined", datasources = datasources)$Global.Mean[, "EstimatedMean"]
  avg_on <- dsBaseClient::ds.mean("df_on_pred", type = "combined", datasources = datasources)$Global.Mean[, "EstimatedMean"]

  if (comparison == "difference") {
    avg_diff <- avg_on - avg_off
  } else if (comparison == "lnriskratio") {
    avg_diff <- log(avg_on / avg_off)
  } else if (comparison == "lnoddsratio") {
    avg_diff <- log((avg_on / (1 - avg_on)) / (avg_off / (1 - avg_off)))
  }
  # avg_diff <- avg_on - avg_off

  return(avg_diff)
}
#' @noRd
meatCL_ds <- function(score_eval,
                      clusters,
                      error_type){

  # calculate the clustered with evals and the corresponding clusters
  n <- NROW(score_eval)
  k <- NCOL(score_eval)
  g <- length(unique(clusters)) # n of clusters

  # HC2/3
  if (error_type == "HC2" | error_type == "HC3"){
    score_eval <- sqrt((g - 1) / g) * score_eval
  }
  ## aggregate within cluster levels
  score_eval_agg <- apply(score_eval, 2L, rowsum, clusters)

  # outer product
  # adj <- 1
  adj <- g / (g - 1)
  rval <- adj * (t(score_eval_agg) %*% score_eval_agg)/n

  # HC0/1
  if (error_type == "HC0"){
    rval <- rval
  }else if (error_type == "HC1"){
    rval <- (n - 1)/(n - k) * rval
  }

  return(rval)
}
#' @noRd
pool_servers <- function(list_servers){
  pool.mat <- matrix(, nrow = 0, ncol = 2)
  for (serv in list_servers){
    pool.mat <- rbind(pool.mat, serv)
  }

  return(pool.mat)
}
#' @noRd
get_range <- function(vec, datasources){
  # get range vals across servers for ecdf

  range_call <- as.symbol(paste0("rangeDS_2(", vec, ")"))
  range_all <- DSI::datashield.aggregate(datasources, range_call)

  range_all <- sapply(range_all, function(x) x)
  out <- c(min(range_all[1,]), max(range_all[2,])) # reintroduce later

  return(out)
}
#' @noRd
get_domain <- function(min, max, len){
  return(seq(min, max, length.out = len))
}
#' @noRd
quantile_from_ecdf <- function(domain, ecdf_vals, probs = seq(0, 1, 0.25)) {
  # sanity checks
  if (length(domain) != length(ecdf_vals))
    stop("domain and ecdf_vals must be same length.")
  if (any(diff(domain) < 0))
    stop("domain must be sorted in increasing order.")
  if (any(ecdf_vals < 0 | ecdf_vals > 1))
    stop("ecdf_vals must lie in [0,1].")
  if (any(probs < 0 | probs > 1))
    stop("probs must lie in [0,1].")

  # Use approx() to invert: treat ecdf_vals as 'x' and domain as 'y'
  out <- stats::approx(
    x    = ecdf_vals,
    y    = domain,
    xout = probs,
    method = "linear",
    ties   = "ordered",
    rule   = 2         # rule=2: if prob<min or >max, returns min(domain)/max(domain)
  )$y

  names(out) <- paste0(probs * 100, "%")
  return(out)
}

# summary functions
#' @noRd
split_summary <- function(unmatched, matched){
  # function for generating summary of unmatched and matched return info

  out_summary <- list()
  cols <- c("Std.Mean.Diff", "Var.Ratio", "eCDF.Max")

  unm <- lapply(unmatched, get_server_summary)
  m <- lapply(matched, get_server_summary)

  out_summary$'Summary of Balance for All Data:' <- unm
  out_summary$'Summary of Balance for Matched Data:' <- m

  # percentage improvement - all except var ratio
  reduction_list <- mapply(
    function(x, y) {out <- 100 * (abs(x[, cols])  - abs(y[, cols])) / abs(x[, cols])
    out[, "Var.Ratio"] <- 100 * (abs(log(x[, "Var.Ratio"]))  - abs(log(y[, "Var.Ratio"]))) / abs(log(x[, "Var.Ratio"]))
    return(out)
    }, unm, m, SIMPLIFY = FALSE)
  out_summary$"Percent Balance Improvement:" <- reduction_list

  return(out_summary)
}
#' @noRd
combine_summary <- function(unmatched, matched){
  # function for generating combined summary of unmatched and matched return info

  out_summary <- list()
  cols <- c("Std.Mean.Diff", "Var.Ratio", "eCDF.Max")

  unm <- get_combined_summary(unmatched)
  m <- get_combined_summary(matched)

  out_summary$'Summary of Balance for All Data:' <- unm
  out_summary$'Summary of Balance for Matched Data:' <- m

  # percentage improvement - all except var ratio
  reduction <- 100 * (abs(unm[, cols])  - abs(m[, cols])) / abs(unm[, cols])
  reduction$Var.Ratio <- 100 * (abs(log(unm[, "Var.Ratio"]))  - abs(log(m[, "Var.Ratio"]))) / abs(log(unm[, "Var.Ratio"]))
  out_summary$"Percent Balance Improvement:" <- reduction

  return(out_summary)
}

# function to assist in summarising
#' @noRd
get_summary <- function(df_sum,
                        df_sumsq,
                        sum_weights){

  df_mean <- df_sum / sum_weights
  df_var <- df_sumsq / (sum_weights-1) - (df_sum^2) / (sum_weights * (sum_weights - 1))

  out <- list(mean = df_mean,
              var = df_var)
  return(out)
}
#' @noRd
get_server_summary <- function(server){
  # function for defining format of summary table for a given server server

  # mean and variance
  control <- get_summary(server$Control$Sum,
                         server$Control$`Sum Squared`,
                         server$Control$sum_weights)
  c_mean <- control$mean
  c_var <- control$var
  treatment <- get_summary(server$Treated$Sum,
                           server$Treated$`Sum Squared`,
                           server$Treated$sum_weights)
  t_mean <- treatment$mean
  t_var <- treatment$var

  # SMD
  smd <- (t_mean - c_mean)
  smd <- smd / sqrt(t_var)

  # var ratio
  var_ratio <- t_var / c_var

  # eCDF
  ecdf_max <- mapply(function(x, y) max(abs(x - y)), server$Control$eCDF, server$Treated$eCDF)

  out <- data.frame("Means.Treated" = t_mean,
                    "Means.Control" = c_mean,
                    "Std.Mean.Diff" = smd,
                    "Var.Ratio" = var_ratio,
                    "eCDF.Max" = ecdf_max
  )

  return(out)
}
#' @noRd
get_combined_summary <- function(df_list){
  # function to aggregate server results

  # get global values from all servers
  global_server <- list()
  global_server$Control$Sum <- global_server$Control$`Sum Squared` <- global_server$Control$sum_weights <- global_server$Control$N <- 0
  global_server$Treated$Sum <- global_server$Treated$`Sum Squared` <- global_server$Treated$sum_weights <- global_server$Treated$N <- 0
  global_server$Control$eCDF <- global_server$Treated$eCDF <- 0

  for (server in df_list){

    # cummulate control vars
    global_server$Control$Sum <- global_server$Control$Sum + server$Control$Sum
    global_server$Control$`Sum Squared` <- global_server$Control$`Sum Squared` + server$Control$`Sum Squared`
    global_server$Control$N <- global_server$Control$N + server$Control$N
    global_server$Control$sum_weights <- global_server$Control$sum_weights + server$Control$sum_weights
    global_server$Control$eCDF <- global_server$Control$eCDF + do.call(cbind.data.frame, server$Control$eCDF)

    global_server$Treated$Sum <- global_server$Treated$Sum + server$Treated$Sum
    global_server$Treated$`Sum Squared` <- global_server$Treated$`Sum Squared` + server$Treated$`Sum Squared`
    global_server$Treated$N <- global_server$Treated$N + server$Treated$N
    global_server$Treated$sum_weights <- global_server$Treated$sum_weights + server$Treated$sum_weights
    global_server$Treated$eCDF <- global_server$Treated$eCDF + do.call(cbind.data.frame, server$Treated$eCDF)
  }

  # # eCDF
  global_server$Control$eCDF <- global_server$Control$eCDF / global_server$Control$sum_weights
  global_server$Treated$eCDF <- global_server$Treated$eCDF / global_server$Treated$sum_weights

  out <- get_server_summary(global_server)

  return(out)
}
#' @noRd
summarise_N <- function(mgroup){
  # function to display sample sizes of the data

  N_control <- 0
  N_treated <- 0
  for (server in mgroup){
    N_control <- N_control + server$Control$N
    N_treated <- N_treated + server$Treated$N
  }

  out <- data.frame(Control = mean(N_control),
                    Treated = mean(N_treated))
  return(out)
}
#' @noRd
reorder_heirarchy <- function(server_list){
  # switch the order of the list from server -> subclass to subclass -> server
  subclass_list <- list()

  # Iterate over each server and its subclasses
  for (server_name in names(server_list)) {
    server <- server_list[[server_name]]

    # Iterate over subclasses within each server
    for (subclass_name in names(server)) {
      subclass_info <- server[[subclass_name]]

      # Check if the subclass already exists in the subclass_list
      if (subclass_name %in% names(subclass_list)) {
        # Append the server information to the existing subclass
        subclass_list[[subclass_name]][[server_name]] <- subclass_info
      } else {
        # Create a new entry for the subclass in the subclass_list
        subclass_list[[subclass_name]][[server_name]] <- subclass_info
      }
    }
  }

  return(subclass_list)
}

#' plotting functions which don't require calls to the servers
#' only matched summary tables required

#' This function creates a Love plot to visualize standardized mean differences before and after matching.
#' It highlights the balance improvement for each covariate in a concise and interpretable way.
#'
#' @param match_obj A list. The output from `ds.match_summary()`,
#' containing both unmatched and matched summary tables.
#'
#' @return A `ggplot` object showing side-by-side standardized mean differences for each covariate
#' across unmatched and matched samples.
#'
#' @details
#' This function extracts balance summary tables before and after matching, reshapes the data for plotting,
#' and generates a Love plot using `ggplot2`. It assumes that the tables include a `Std.Mean.Diff` column
#' and named rows corresponding to covariates.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item None. This is a purely client-side visualization function.
#' }
#'
#' @author Roy Gusinow
#'
#' @export
love.plot <- function(match_obj){
  # unmatched, matched are tables from the summary table

  unmatched <- match_obj$`Summary of Balance for All Data:`
  matched <- match_obj$`Summary of Balance for Matched Data:`


  unmatched$type <- "All Data"
  unmatched$var <- rownames(unmatched)
  matched$type <- "Matched Data"
  matched$var <- rownames(matched)
  df <- rbind(unmatched, matched)


  text.size = 14
  frame <- ggplot(df, aes(x = !!sym("Std.Mean.Diff"), y = !!sym("var"), color = type)) +
    geom_point(size = 5) +
    geom_vline(xintercept = 0) +

    theme(panel.background = element_rect(fill = "white"),

          axis.text.x = element_text(color = "black", size = text.size - 1),
          axis.text.y = element_text(color = "black", size = text.size - 1),
          axis.title= element_text(size = text.size + 1, face="italic"),
          plot.title = element_text(size=text.size + 3, face = "bold"),

          panel.border = element_rect(fill = NA, color = "black"),
          plot.background = element_blank(),

          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title=element_text(size=text.size - 2),
          legend.text=element_text(size=text.size - 3)) +


    scale_fill_discrete(name = "Sample") +
    ggtitle("Covariate Balance") +
    xlab("Standardized Mean Differences") +
    ylab("Variables")

  return(frame)
}

#' This function creates jittered dot plots of propensity score distributions
#' for treated and control units, before and after matching. It provides a visual
#' comparison of how the distributions align across groups.
#'
#' @param match A list. The result of a `ds.match_summary()` call.
#'
#' @return A faceted `ggplot` grid combining four panels:
#' \itemize{
#'   \item Unmatched Treated Units
#'   \item Matched Treated Units
#'   \item Matched Control Units
#'   \item Unmatched Control Units
#' }
#' Each panel shows a jittered horizontal layout of propensity scores.
#'
#' @details
#' This function visualizes the overlap of propensity scores between treated and control
#' groups, before and after matching. Treated and control units are plotted separately,
#' with jitter applied to reduce overplotting. It helps assess the quality of matching
#' in terms of covariate overlap.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item None. This is a client-side visualization function using only local match output.
#' }
#'
#' @author Roy Gusinow
#' @export
jitter.plot <- function(match){
  # requires the list of propensity scores currently
  # neds checking

  unmatched <- match$prop_pool
  matched <- match$pooled.match

  x <- unmatched[, "distance"]
  y = rep(1, length(x))
  y_jitter <- jitter(y)

  get_plot <- function(df, treat = 1, title, keep_x = F){

    p <- ggplot(df[df$treat == treat,], aes(x = !!sym("distance"), y = jitter(rep(1, length(distance))))) +
      geom_point(size = 5, alpha = 1/10) +
      xlim(0,1) +
      ylim(0.95,1.05) +

      theme_minimal() +
      ggtitle(title) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.title.y = element_blank(), axis.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5))

    if(keep_x){
      p <- p + theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size=10)) + xlab("Propensity Score")
    }
    return(p)
  }

  p1 <- get_plot(unmatched, treat = 1, "Unmatched Treated Units")
  p2 <- get_plot(matched, treat = 1, "Matched Treated Units")
  p3 <- get_plot(matched, treat = 0, "Matched Control Units")
  p4 <- get_plot(unmatched, treat = 0, "Unmatched Control Units", keep_x = T)

  plot_grob <- rbind(
    ggplotGrob(p2),
    ggplotGrob(p3),
    ggplotGrob(p1),
    ggplotGrob(p4),
    size = "last"
  )

  grid.newpage()
  grid.draw(plot_grob)

  return(plot_grob)

}

#' @noRd
plot_single_qq <- function(df){
  # plot a single qq
  # df should contain 2 columns
  # qq_untreated and qq_treated

  rr <- range(c(df$qq_control, df$qq_treated))
  plt <-
    ggplot(df, aes(x = !!sym("qq_control"), y = !!sym("qq_treated"))) +
    geom_abline(intercept = 0, slope = 1, color="red") +
    geom_point() +
    xlim(rr) + ylim(rr) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


  return(plt)
}
#' @noRd
get_qq_plot_per_var <- function(var, datasources){
  qq_control <- dsBaseClient::ds.quantileMean(x = paste0("control$", var),
                                type = "combine", datasources = datasources)[1:7]
  qq_treated <- dsBaseClient::ds.quantileMean(x = paste0("treated$", var),
                                type = "combine", datasources = datasources)[1:7]

  plt <- plot_single_qq(data.frame(qq_control = qq_control, qq_treated = qq_treated)) + ylab(var)

  return(plt)
}

# eCDF
#' @noRd
get_eCDF_plot_per_var <- function(var, len, datasources){
  # get the grouped ecdf plot

  ecdf_control <- ds.eCDF(object = paste0("control$", var),
                          type = "combine",
                          len = len,
                          datasources = datasources)
  ecdf_treated <- ds.eCDF(object = paste0("treated$", var),
                          type = "combine",
                          len = len,
                          datasources = datasources)

  control_df <- data.frame(ecdf_control)
  treated_df <- data.frame(ecdf_treated)
  control_df$type <- "Control"
  treated_df$type <- "Treated"
  plt <- plot_single_eCDF(rbind(control_df, treated_df)) + ylab(var)

  return(plt)
}
#' @noRd
plot_single_eCDF <- function(df){
  # plot a single eCDF
  # df should contain 2 columns
  # ecdf_untreated and ecdf_treated

  rr_domain <- range(df[,"domain"])
  rr_ecdf <- range(df[,"ecdf"])
  plt <-
    ggplot(df, aes(x = !!sym("domain"), y = !!sym("ecdf"), color = type)) +
    geom_step() +
    xlim(rr_domain) + ylim(rr_ecdf) +

    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    geom_hline(yintercept=1, linetype="dashed", color = "black") +

    scale_color_manual(values=c("#11468F", "#DA1212")) +

    theme_minimal() +
    theme(axis.title.x=element_blank(),

          legend.position="none")


  return(plt)
}


# misc
#' @noRd
ds.delete_col <- function(
  data,
  col_name,
  datasources = NULL
){
  # delete specified column from ds frame if there

  assign_ONES(paste0(data, "$", col_name), datasources = datasources)

  # 2. Retrieve column names and find the index of "distance"
  colnames_list <- dsBaseClient::ds.colnames(
    x           = data,
    datasources = datasources
  )

  server1_names <- colnames_list$server1
  idx <- which(server1_names == col_name)

  if (length(idx) != 0){
    # 3. Subset DST to drop that column
    dsBaseClient::ds.dataFrameSubset(
      df.name         = data,
      V1.name         = "ONES",
      V2.name         = "ONES",
      Boolean.operator= "==",
      keep.cols       = NULL,
      rm.cols         = idx,
      newobj          = data,
      datasources     = datasources
    )
  }

}

#' @noRd
remove_weights <- function(formula) {
  # remove weights from formula object
  fstr <- paste(deparse(formula), collapse = "")
  fstr <- gsub("\\+?\\s*weights\\([^)]*\\)", "", fstr)
  fstr <- gsub("\\s+", " ", fstr)
  fstr <- trimws(fstr)
  fstr <- sub('^"(.*)"$', '\\1', fstr) # remove leading/trailing quotes if present
  formula_clean <- as.formula(fstr, env = parent.frame())
  return(formula_clean)
}


#' @noRd
encode_coef_names <- function(coef_vec) {
  names(coef_vec) <- gsub(":", "___", names(coef_vec), fixed = TRUE)
  coef_vec
}

# strip_weights_from_formula <- function(formula_obj) {
#   # Convert to character for easy manipulation
#   # formula_str <- deparse(formula_obj)
#   # Collapse multi-line formula if necessary
#   formula_str <- paste(formula_obj, collapse = " ")
#   # Remove any '+ weights(...)' term (with or without surrounding spaces)
#   formula_str <- sub("\\+\\s*weights\\s*\\([^\\)]*\\)\\s*", "", formula_str)
#   # Convert back to formula
#   return (formula_str)
# }

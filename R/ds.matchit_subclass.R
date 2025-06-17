#' This function performs subclassification-based propensity score matching in a federated setting.
#' It estimates propensity scores via a logistic regression, computes quantiles of the scores based
#' on the specified estimand (ATE, ATT, or ATC), and assigns subclass labels accordingly.
#'
#' @param formula A formula object. Specifies the GLM used to estimate the propensity score. The treatment variable must appear on the left-hand side.
#' @param data A character string. The name of the data frame on the server-side.
#' @param subclass An integer. Number of subclasses to create by slicing the propensity score distribution. Default is `6`.
#' @param estimand A character string. The causal estimand of interest: `"ATE"`, `"ATT"`, or `"ATC"`. Default is `"ATE"`.
#' @param newobj A character string. The name of the new object to create on each server containing the matched dataset. Default is `"matched_df_list"`.
#' @param len An integer. The number of bins used in the empirical CDF for quantile estimation. Default is `429`.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function uses `datashield.connections_find()` to detect active connections.
#'
#' @return A list containing:
#' \describe{
#'   \item{Summary}{Balance summaries before and after subclassification matching.}
#'   \item{Sample Sizes}{A table showing treated/control counts within each subclass.}
#'   \item{quantiles}{The estimated propensity score quantile cut points used for subclassification.}
#'   \item{subclass}{A numeric vector indexing subclass IDs.}
#' }
#'
#' @details
#' This function runs a logistic GLM to estimate propensity scores, uses pooled ECDFs to determine subclass thresholds,
#' and then assigns subclass labels accordingly using `matchit_subclassDS`. It also computes and applies
#' subclass weights depending on the estimand.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item `matchit_subclassDS`
#'   \item `assign_weights_subclassDS`
#'   \item `rangeDS`
#' }
#'
#' @author Roy Gusinow

#' @export
ds.matchit_subclass <- function(formula,
                                data,
                                subclass = 6,
                                estimand = "ATE",
                                newobj = "matched_df",
                                len = 429,
                                datasources = NULL
){

  # look for DS datasources
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }

  # if the argument 'data' is set, check that the data frame is defined (i.e. exists) on the server site
  if((is.null(data))){
    stop("Please provide a valid dataframe on the serverside")
  }

  treatment = all.vars(formula)[1] # get treatment variable

  # generate propensity on the server
  fed_glm <- dsBaseClient::ds.glm(formula = formula,
                    data = data,
                    family = "binomial",
                    datasources = datasources)

  # gen prop scores on the server side - using fed coeffs
  coef_vec <- fed_glm$coefficients[, "Estimate"]
  names(coef_vec)[names(coef_vec) == "(Intercept)"] <- "Intercept"

  ds.genProp(formula = formula, coefficients = coef_vec, data = data, newobj = "distance", link = "logit", datasources = datasources)
  dsBaseClient::ds.cbind(x = c(data, "distance"),
           newobj = data,
           datasources = datasources)

  # subclass start
  # get quantiles from eCDF
  if (estimand == "ATT"){
    assign_ONES(paste0(data, "$", treatment), datasources = datasources)
    dsBaseClient::ds.dataFrameSubset(df.name = data,
                       V1.name = paste0(data, "$", treatment),
                       V2.name = "ONES",
                       Boolean.operator = "==",
                       keep.cols = NULL, #All columns are included in the new subset
                       rm.cols = NULL, #All columns are included in the new subset
                       keep.NAs = FALSE, #All rows with NAs are removed
                       newobj = "temp",
                       datasources = datasources,#only the first server is used ("study1")
                       notify.of.progress = FALSE)
    ecdf_x <- ds.eCDF("temp$distance",
                      type = "combine",
                      len = len,
                      datasources = datasources
    )
  } else if (estimand == "ATC"){
    assign_ZEROS(paste0(data, "$", treatment), datasources = datasources)
    dsBaseClient::ds.dataFrameSubset(df.name = data,
                       V1.name = paste0(data, "$", treatment),
                       V2.name = "ZEROS",
                       Boolean.operator = "==",
                       keep.cols = NULL, #All columns are included in the new subset
                       rm.cols = NULL, #All columns are included in the new subset
                       keep.NAs = FALSE, #All rows with NAs are removed
                       newobj = "temp",
                       datasources = datasources,#only the first server is used ("study1")
                       notify.of.progress = FALSE)
    ecdf_x <- ds.eCDF("temp$distance",
                      type = "combine",
                      len = len,
                      datasources = datasources
    )
  } else if (estimand == "ATE"){
    ecdf_x <- ds.eCDF(paste0(data, "$distance"),
                      type = "combine",
                      len = len,
                      datasources = datasources
    )
  }

  seq_probs <- seq(0, 1, length.out = round(subclass) + 1)
  quantiles <- quantile_from_ecdf(ecdf_x$domain, ecdf_x$ecdf, probs = seq_probs)
  names(quantiles) <- NULL

  cally2 <- call("matchit_subclassDS",
                 data = data,
                 distance = paste0(data, "$distance"),
                 quantiles = quantiles)
  result <- DSI::datashield.assign.expr(datasources, newobj, cally2)

  # get weights
  fed_summary.combined <- ds.match_summary_subclass(unmatched_obj = data,
                                                    matched_obj = newobj,
                                                    # type = "combined",
                                                    # bin_num = 429,
                                                    treatment = treatment,
                                                    datasources = datasources)
  ss <- fed_summary.combined$`Sample Sizes`
  if (estimand == "ATT"){
    treat_over_control <- ss[,2] / ss[,1]
    names(treat_over_control) <- paste0("T", 1:length(treat_over_control))
    weights_list <- treat_over_control
  } else if (estimand == "ATC") {
    control_over_treat <- ss[,1] / ss[,2] # atc
    names(control_over_treat) <- paste0("C", 1:length(control_over_treat))
    weights_list <- control_over_treat
  } else if (estimand == "ATE") {
    total_treat_ratio <- ss[,3] / ss[,2]
    total_control_ratio <- ss[,3] / ss[,1]

    names(total_treat_ratio) <- paste0("T", 1:length(total_treat_ratio))
    names(total_control_ratio) <- paste0("C", 1:length(total_control_ratio))
    weights_list <- c(total_treat_ratio, total_control_ratio)
  }

  # delete weights col
  ds.delete_col(
    data = newobj,
    col_name = "weights",
    datasources = datasources
  )

  # assign weights
  cally3 <- call("assign_weights_subclassDS",
                 data = newobj,
                 weights_list = weights_list,
                 treatment = treatment,
                 estimand = estimand)
  assign_weights <- DSI::datashield.assign.expr(datasources, "weights", cally3)
  dsBaseClient::ds.cbind(x = c(newobj, "weights"),
           newobj = newobj,
           datasources = datasources)

  fed_summary.combined$quantiles <- quantiles
  fed_summary.combined$subclass <- 1:sum(ss[,3])

  return(fed_summary.combined)
}

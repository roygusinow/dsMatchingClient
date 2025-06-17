#' This function generates subclass-specific covariate balance summaries after matching.
#' It compares covariate distributions between treated and control groups for each subclass,
#' optionally reporting percent improvements relative to the unmatched sample.
#'
#' @param unmatched_obj A character string. The name of the unmatched (original) server-side data frame.
#' @param matched_obj A character string. The name of the matched server-side data frame that contains subclass assignments.
#' @param detailed A logical value. Whether to return percent improvements in balance measures (standardized mean difference, variance ratio, ECDF maximum) for each subclass. Default is `FALSE`.
#' @param bin_num An integer. Number of bins to use in ECDF binning. Default is `40`.
#' @param treatment A character string. The name of the treatment variable in the server-side data.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function uses `datashield.connections_find()` to detect active connections.
#'
#' @return A named list containing:
#' \describe{
#'   \item{Summary}{A nested list with summaries for unmatched data and each matched subclass.}
#'   \item{Percent Balance Improvement (optional)}{If `detailed = TRUE`, includes percent change in balance metrics relative to unmatched data.}
#'   \item{Sample Sizes}{A table showing the number of treated and control units in each subclass.}
#' }
#'
#' @details
#' This function allows for fine-grained diagnostics of covariate balance after subclassification-based
#' matching in a federated DataSHIELD setting. It bins covariates using pooled ranges across studies
#' and excludes non-numeric or categorical variables.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item `match_summaryDS`
#'   \item `match_summary_subclassDS`
#'   \item `rangeDS`
#' }
#'
#' @author Roy Gusinow
#'
#' @export
ds.match_summary_subclass <- function(unmatched_obj,
                                      matched_obj,
                                      detailed = F,
                                      bin_num = 40,
                                      treatment,
                                      datasources = NULL

){

  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  cnames <- dsBaseClient::ds.colnames(unmatched_obj, datasources = datasources)[[1]]

  # remove factor char cols + others
  class_col_mask <- sapply(X = as.list(cnames),
                           FUN = function(vari, df_name){
                             class_name <- dsBaseClient::ds.class(paste0(df_name, "$", vari), datasources = datasources)[[1]]
                             return(class_name == "factor" | class_name == "character")},
                           df_name = unmatched_obj, USE.NAMES = TRUE, simplify = T)
  factor_cnames <- cnames[class_col_mask]
  cnames <- cnames[! cnames %in% c( "distance", factor_cnames, treatment)] #"weights",

  # get range vals across servers for ecdf
  cnames <- c("distance", cnames)
  range_names <- paste0(unmatched_obj, "$", cnames)
  range_frame <- sapply(X = as.list(range_names), FUN = get_range, datasources = datasources, USE.NAMES = TRUE, simplify =T)

  # convert to char format
  range_frame.char <- paste0(as.character(rbind(cnames, range_frame)), collapse = ",")

  # first call for unmatched
  cally <- call("match_summaryDS", unmatched_obj, bin_num, treatment, range_frame.char, weights = FALSE)
  unmatched <- DSI::datashield.aggregate(datasources, cally)

  # first call for matched
  cally2 <- call("match_summary_subclassDS", matched_obj, bin_num, treatment, range_frame.char, weights = TRUE)
  matched_list <- DSI::datashield.aggregate(datasources, cally2)

  # generate summary based on type
  reordered_match <- reorder_heirarchy(matched_list)
  sub_list <- list()
  sub_no <- 1
  for (subclass in reordered_match){
    m_sum <- get_combined_summary(subclass)
    sub_list[[paste0("Subclass ", sub_no)]] <- m_sum
    sub_no <- sub_no + 1
  }

  out_summary <- list()
  unm <- get_combined_summary(unmatched)
  out_summary$'Summary'$'Summary of Balance for All Data:' <- unm
  out_summary$'Summary'$'Summary of Balance for Matched Data:' <- sub_list

  # add detailed table
  if (detailed){
    cols <- c("Std.Mean.Diff", "Var.Ratio", "eCDF.Max")
    reduction_list <- mapply(
      function(x, y) {out <- 100 * (abs(x[, cols])  - abs(y[, cols])) / abs(x[, cols])
      out[, "Var.Ratio"] <- 100 * (abs(log(x[, "Var.Ratio"]))  - abs(log(y[, "Var.Ratio"]))) / abs(log(x[, "Var.Ratio"]))
      return(out)
      }, rep(list(unm), length(sub_list)), sub_list, SIMPLIFY = FALSE, USE.NAMES = T)
    names(reduction_list) <- paste0(names(sub_list), " Improvement")
    out_summary$"Percent Balance Improvement:" <- reduction_list
  }

  ls_N <- lapply(reordered_match, summarise_N)
  s_tot <- do.call("rbind", ls_N)
  s_tot$"Total" <- s_tot$Control + s_tot$Treated
  rownames(s_tot) <- paste0("Subclass ", names(reordered_match))
  out_summary$'Sample Sizes' <- s_tot

  return(out_summary)
}














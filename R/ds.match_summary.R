#' This function generates a covariate balance summary table before and after matching.
#' It compares covariate distributions between treated and control groups using empirical
#' cumulative distribution functions (ECDFs), with optional binning and summary type.
#'
#' @param unmatched_obj A character string. The name of the unmatched (original) server-side data frame.
#' @param matched_obj A character string. The name of the matched server-side data frame produced by matching (e.g., `ds.matchit`).
#' @param type A character string. The type of summary to return: `"split"` for study-wise results or `"combined"` for pooled summaries. Default is `"combine"`.
#' @param bin_num An integer. Number of bins to use in ECDF binning. Default is `40`.
#' @param treatment A character string. The name of the treatment variable in the server-side data.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function uses `datashield.connections_find()` to detect active connections.
#'
#' @return A named list containing covariate balance summaries. Includes pre- and post-matching comparisons
#' for each covariate and a sample size table with matched and unmatched counts.
#'
#' @details
#' This function is used to evaluate balance in observed covariates between treated and control units
#' before and after matching in federated DataSHIELD settings. Factor and character variables are excluded.
#' Balance is assessed via ECDF-based summaries using binning across a range pooled across all studies.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item `match_summaryDS`
#'   \item `rangeDS`
#' }
#'
#' @author Roy Gusinow
#'
#' @export
ds.match_summary <- function(unmatched_obj,
                             matched_obj,
                             type = "combine",
                             bin_num = 40,
                             treatment = NULL,
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
  print(cally)
  unmatched <- DSI::datashield.aggregate(datasources, cally)

  # first call for matched
  cally2 <- call("match_summaryDS", matched_obj, bin_num, treatment, range_frame.char, weights = TRUE)
  matched <- DSI::datashield.aggregate(datasources, cally2)

  # generate summary based on type
  if(type == "split"){
    out_summary <- split_summary(unmatched, matched)

  }else if(type == 'combined'){
    out_summary <- combine_summary(unmatched, matched)
  }

  # add sample size info
  s1 <- summarise_N(unmatched)
  s2 <- summarise_N(matched)
  out_summary$'Sample Sizes' <- rbind(s1, s2)
  rownames(out_summary$'Sample Sizes') <- c("All", "Matched")

  return(out_summary)
}

#' Compute Empirical Cumulative Distribution Function (eCDF) for a Variable in DataSHIELD
#'
#' This function computes the empirical cumulative distribution function (eCDF) for a numeric variable
#' using a server-side aggregate function `eCDFDS`. It returns either split (study-specific) or combined
#' (pooled) normalized results.
#'
#' @param object A character string. Name of the variable (e.g., `"DST$age"`) on the server-side to compute the eCDF for.
#' @param type A character string. One of `"split"` (default) to return normalized eCDFs per study, or `"combine"` to return a pooled eCDF.
#' @param len An integer. Number of points in the domain over which the eCDF is evaluated. Default is 40.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function attempts to find active connections via `datashield.connections_find()`.
#'
#' @return A named list with:
#' \describe{
#'   \item{domain}{Numeric vector of domain points for the eCDF.}
#'   \item{ecdf}{Numeric vector of the corresponding eCDF values.}
#' }
#'
#' @details
#' This function is used internally to produce smooth eCDF plots for matched vs unmatched comparisons.
#'
#' #' **Server-side functions called**:
#' #' \itemize{
#' #'   \item `eCDFDS`
#' #' }
#'
#'  @author Roy Gusinow
#'
#' @export
ds.eCDF <- function(object,
                    type = "combine",
                    len = 40,
                    datasources = NULL

){

  # look for DS connections
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }
  if (!type %in% c("split", "combine")) {
    stop("Parameter 'type' must be either 'split' or 'combine'.")
  }
  if (!is.numeric(len) || length(len) != 1 || len <= 0) {
    stop("Parameter 'len' must be a single positive numeric value.")
  }

  object_range <- get_range(object, datasources = datasources)
  domain <- get_domain(min = object_range[1],
                       max = object_range[2],
                       len = len)

  cally <- call("eCDFDS",
                object = object,
                min = object_range[1],
                max = object_range[2],
                len = len,
                weights = FALSE) # weights currently useless
  rval <- DSI::datashield.aggregate(datasources, cally)

  if (type == "split"){
    max_list <- lapply(rval, max)
    eCDF <- mapply(function(x, y){x / y}, rval, max_list)
  }else if (type == "combine"){
    eCDF <- Reduce("+", rval)
    eCDF <- eCDF / max(eCDF)
  }

  # at first element for x axis
  domain <- c(domain[1] - domain[1] * 0.001, domain)
  eCDF <- c(0, eCDF)
  out_list <- list(domain = domain,
                   ecdf = eCDF)

  return(out_list)
}

#' Assign Hat Matrix Diagonal to Server-side Object
#'
#' Computes and assigns the diagonal of the hat matrix (leverages) from a GLM model to each study server.
#' These diagonals are useful for influence diagnostics and robust error estimation.
#'
#' @param fit A list. The GLM output object from `ds.glm`, which must contain at least a `formula` and sample size (`Ntotal`).
#' @param data A character string. The name of the server-side data frame object.
#' @param newobj_name A character string. The name of the new object to be created on the server. Default is `"h_diag"`.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, will attempt to find active connections via `datashield.connections_find()`.
#'
#' @return No return value
#'
#' @details
#' In federated settings using DataSHIELD, hat matrix diagonals cannot be retrieved directly but can be reconstructed
#' using global covariate means and the global covariance matrix. This function coordinates that process and assigns
#' the result server-side for later use in robust variance calculations.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item `hat_diagDS`
#' }
#'
#'
#' @author Roy Gusinow
#'
#' @export
ds.hat_diag <- function(fit,
                        data,
                        newobj_name = "h_diag",
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

  if((is.null(fit)) | !is.list(fit)){
    stop("Please provide a valid GLM output list from ds.glm")
  }

  if((is.null(data))){
    stop("Please provide a valid dataframe on the serverside")
  }

  # Ensure fit$formula is a formula; convert if not
  if (!inherits(fit$formula, "formula")) {
    fit_formula <- stats::as.formula(fit$formula)
  }

  # get global components
  # cov mat
  global_cov_mat <- dsBaseClient::ds.cov(x = data,
                           type = "combine",
                           datasources = datasources)$`Variance-Covariance Matrix`
  # avg of vars
  global_avg_vec <- c()
  cnames <- dsBaseClient::ds.colnames(data, datasources = datasources)[[1]]
  for (var_name in cnames){
    var_mean <- dsBaseClient::ds.mean(x = paste0(data, "$", var_name), type = "combine", datasources = datasources)$Global.Mean[, "EstimatedMean"]
    global_avg_vec <- c(global_avg_vec, var_mean)
  }
  names(global_avg_vec) <- cnames

  cally <- call("hat_diagDS",
                form = fit_formula, data,
                global_cov_mat = c(global_cov_mat),
                global_avg_vec,
                global_n = fit$Ntotal)
  DSI::datashield.assign.expr(datasources, newobj_name, cally)

  return(NULL)
}

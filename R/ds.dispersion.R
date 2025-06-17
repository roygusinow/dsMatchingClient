#' Compute Dispersion Estimate for GLM (Used in Sandwich Estimator)
#'
#' This function computes the dispersion (variance inflation) term for a GLM model
#' fitted via DataSHIELD. It is used internally for robust standard error estimation,
#' especially with sandwich (HC-type) estimators.
#'
#' @param fit A list. Output from `ds.glm`, must include a `formula` and `coefficients`.
#' @param data A character string. Name of the data frame stored on the server side.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, attempts to find active connections using `datashield.connections_find()`.
#'
#' @return A numeric value. The estimated dispersion across studies (pooled).
#'
#' @details
#' This function generates predicted values using the GLM link function and then computes
#' a dispersion estimate by comparing residual and fitted values aggregated across studies.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item `dispersionDS`
#' }
#'
#'  @author Roy Gusinow
#'
#' @export
ds.dispersion <- function(fit,
                          data,
                          datasources = NULL

){
  # look for DS connections
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  # Input checks
  if (!(is.list(datasources) &&
        all(vapply(datasources, function(d) methods::is(d, "DSConnection"), logical(1))))) {
    stop("The 'datasources' argument must be a list of DSConnection-class objects.", call. = FALSE)
  }

  if (!is.list(fit) || is.null(fit$formula)) {
    stop("The 'fit' argument must be a list (from ds.glm) containing at least a 'formula' component.")
  }

  if (!is.character(data) || length(data) != 1) {
    stop("The 'data' argument must be a character string referring to a server-side data object.")
  }

  # Ensure fit$formula is a formula; convert if not
  if (!inherits(fit$formula, "formula")) {
    formula <- stats::as.formula(fit$formula)
  }

  # prediction vector
  coef_vec <- fit$coefficients[,1]
  names(coef_vec)[names(coef_vec) == "(Intercept)"] <- "Intercept"

  ds.genProp(formula = formula, coefficients = coef_vec, data = data, newobj = "pred_col", link = fit$family$link, datasources = datasources)

  # dispersion call
  cally <- call("dispersionDS", formula, data, "pred_col", fit$family$family)
  dispersion_list <- DSI::datashield.aggregate(datasources, cally)

  # weights call
  mean_list <- dsBaseClient::ds.mean(paste0(data, "$weights"), type = "split", datasources = datasources)
  sum_list <- mapply(function(x, y) x * y, mean_list$Mean.by.Study[, "EstimatedMean"], mean_list$Mean.by.Study[, "Ntotal"])

  dispersion <- Reduce("+", dispersion_list) / Reduce("+", sum_list)

  return(dispersion)

}

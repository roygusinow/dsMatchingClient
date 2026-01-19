#' This function evaluates the estimating (score) function for each observation in a GLM model.
#' In the simplest case, this corresponds to the product of the design matrix and residuals, i.e.,
#' the derivative of the log-likelihood with respect to parameters. It is primarily used for
#' heteroskedasticity-consistent (HC) variance estimation.
#'
#' @param fit A list. Output from `ds.glm`, containing at least `formula`, `coefficients`, and `family`.
#' @param data A character string. The name of the data frame on the server-side.
#' @param error_type A character string specifying the HC-type estimator to use. One of `"HC0"`, `"HC1"`, `"HC2"`, or `"HC3"`. Default is `"HC0"`.
#' @param h_diag A character string. The name of the object on the server that contains the diagonals of the hat matrix. Typically created with `ds.hat_diag()`.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function uses `datashield.connections_find()` to detect active connections.
#'
#' @return A list of matrices (one per study) containing scaled estimating functions for each observation and model parameter.
#'
#' @details
#' In DataSHIELD it is not possible to see the data on the servers of the collaborating studies.
#' It is only possible to get summaries of objects stored on the server-side. This function computes
#' the estimating function (score function) using predicted values and hat matrix diagonals, for use in
#' robust variance estimation.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item `estfunDS`
#' }
#'
#'  @author Roy Gusinow
#'
#' @export
ds.estfun <- function(fit,
                      data,
                      error_type = "HC0",
                      h_diag = "h_diag",
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

  if (!error_type %in% paste0("HC", 0:3)) {
    stop("Please choose either HC0, HC1, HC2, or HC3 for the error specification.")
  }

  if((is.null(h_diag))){
    stop("Please provide an object name on the server side for the hat diagonals. Consider using ds.hat_diag")
  }

  fit_formula <- stats::as.formula(fit$formula)

  # prediction vector
  coef_vec <- fit$coefficients[,1]
  names(coef_vec)[names(coef_vec) == "(Intercept)"] <- "Intercept"
  ds.genProp(formula = fit_formula, coefficients = coef_vec, data = data, newobj = "pred_col", link = fit$family$link, datasources = datasources)

  # dispersion
  disp <- ds.dispersion(fit,
                        data,
                        datasources = datasources)

  # estfun call
  formula_clean <- remove_weights(fit_formula)

  cally <- call("estfunDS", formula_clean, data, "pred_col", error_type, h_diag)
  rval_unscaled_list <- DSI::datashield.aggregate(datasources, cally)

  rval_scaled_list <- mapply(function(x, y) x / y,
                             rval_unscaled_list,
                             rep(disp, length(rval_unscaled_list)))

  return(rval_scaled_list)
}

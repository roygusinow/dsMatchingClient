#' Generate Propensity Scores on the Server Side
#'
#' @param formula A formula object. The regression formula used to compute the scores.
#' @param coefficients A numeric vector. Coefficients corresponding to the variables in the formula.
#' @param data A character string. The name of the data frame on the server side to operate on.
#' @param newobj A character string. The name of the new object to be created on the server side. Default is `"propscore"`.
#' @param link A character string. The link function used in the regression (e.g., `"logit"`, `"identity"`). Default is `"identity"`.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function uses `datashield.connections_find()` to detect active connections.
#'
#' @return No return value. A new object is created on the server side.
#'
#' @details
#' In DataSHIELD it is not possible to see the data on the servers of the collaborating studies. It is only possible to get summaries of objects stored on the server-side. This function creates a new server-side object containing the predicted values from a linear predictor using a specified regression formula and coefficients.
#'
#' Server function called: genPropDS
#'
#' @author Roy Gusinow
#'
#' @export
ds.genProp <- function(formula,
                       coefficients,
                       data,
                       newobj = "propscore",
                       link = "identity",
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

  check_formula(formula)

  # coefficient check
  if(!inherits(coefficients, "numeric")){
    stop("Please use a numeric vector for coefficients", call.=FALSE)
  }

  # if the argument 'data' is set, check that the data frame is defined (i.e. exists) on the server site
  if((is.null(data))){
    stop("Please provide a valid dataframe on the serverside")
  }

  # weight  removal
  formula_clean <- remove_weights(formula)
  coefficients <- encode_coef_names(coefficients)

  # first call
  cally <- call("genPropDS", formula_clean, coefficients, data, link)
  result <- DSI::datashield.assign.expr(datasources, newobj, cally)

  return()
}

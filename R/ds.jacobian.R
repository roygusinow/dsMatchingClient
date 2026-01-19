#' This function computes the Jacobian of the average treatment effect (ATE, ATT, or ATC)
#' with respect to the model parameters in a GLM. It is used in the computation of heteroskedasticity-consistent
#' variance estimators in federated settings.
#'
#' @param fit A list. Output from a `ds.glm` model. Must contain required elements like the formula and coefficients.
#' @param data A character string. The name of the data frame stored on the server side.
#' @param treatment A character string. The name of the treatment variable within `data`.
#' @param estimand A character string specifying the estimand to compute. One of `"ATE"`, `"ATT"`, or `"ATC"`. Default is `"ATE"`.
#' @param eps A numeric value. The finite difference step size for approximating the Jacobian. Default is `1e-7`.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function uses `datashield.connections_find()` to detect active connections.
#'
#' @return A numeric matrix representing the Jacobian of the treatment effect with respect to the model parameters.
#'
#' @details
#' In federated analyses, the Jacobian is approximated using finite differences by perturbing
#' each coefficient in the GLM and recomputing the treatment effect.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item None
#' }
#'
#' @author Roy Gusinow
#'
#' @export
ds.jacobian <- function(fit,
                        data,
                        treatment,
                        estimand = "ATE",
                        comparison = "difference",
                        eps = 1e-7,
                        datasources = NULL

){

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

  if (!is.character(treatment) || length(treatment) != 1) {
    stop("The 'treatment' argument must be a character string indicating a variable in the data.")
  }

  if (!estimand %in% c("ATE", "ATT", "ATC")) {
    stop("The 'estimand' must be one of: 'ATE', 'ATT', or 'ATC'.")
  }

  if (!comparison %in% c("difference", "lnriskratio", "lnoddsratio")) {
    stop("The 'comparison' must be one of: 'difference' or 'ratio'.")
  }

  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) {
    stop("The 'eps' argument must be a single positive numeric value.")
  }

  # get ATE
  clean_avg <- avg_compute_estimate(fit,
                                    data,
                                    treatment,
                                    estimand,
                                    comparison,
                                    datasources)

  # construct J
  J <- c()

  for (i in 1:length(fit$coefficients[, "Estimate"])){
    perturbed_fit <- fit

    fit_coef_perturbed <- perturbed_fit$coefficients[, "Estimate"]
    fit_coef_perturbed[i] <- fit_coef_perturbed[i] + eps
    perturbed_fit$coefficients[, "Estimate"] <- fit_coef_perturbed

    perturbed_avg <- avg_compute_estimate(perturbed_fit,
                                          data,
                                          treatment,
                                          estimand,
                                          comparison,
                                          datasources)

    J <- c(J, (clean_avg - perturbed_avg) / eps)
  }
  J <- as.matrix(J)

  return(J)

}

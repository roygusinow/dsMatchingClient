#' Compute Average Treatment Effects (ATE, ATT, ATC) with Optional Standard Errors
#'
# This function computes the average treatment effect (ATE), average treatment effect on the treated (ATT),
# or average treatment effect on the controls (ATC) from a fitted DataSHIELD GLM object. It can optionally return
# heteroskedasticity-consistent standard errors using sandwich variance estimation and cluster-robust correction.
#'
#' @param fit A list. Output from a `ds.glm` model. Must contain required elements like the formula and coefficients.
#' @param data A character string. The name of the data frame stored on the server side.
#' @param treatment A character string. The name of the treatment variable within `data`.
#' @param estimand A character string specifying the estimand to compute. One of `"ATE"`, `"ATT"`, or `"ATC"`. Default is `"ATE"`.
#' @param comparison A character string indicating the type of comparison. One of `"difference"` (default), `"lnriskratio"`, or `"lnoddsratio"`.
#' @param error_type A character string indicating the HC-type estimator to correct bias (e.g., `"HC0"`, `"HC1"`, etc.).
#' If `NULL`, standard errors are not returned.
#' @param clusters A factor vector indicating cluster membership (e.g., subclass assignments). Typically obtained via `ds.matchit`.
#' @param eps A numeric value for the finite difference step size used in Jacobian approximation. Default is `1e-7`.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, will attempt to find active connections via `datashield.connections_find()`.
#'
#' @return A named list with at least:
#' \describe{
#'   \item{estimate}{The computed average treatment effect.}
#'   \item{se}{(Optional) The standard error of the estimate, if `error_type` is provided.}
#' }
#'
#' @details
#' This function wraps internal client-side calls to compute treatment effects and, optionally, robust standard errors.
#' It is part of the `dsMatchingClient` workflow and expects GLM models fitted using `ds.glm`.
#'
#' After the propensity score matching procedure has been conducted, there exist
#' \eqn{H} subclassification strata (or matched groups) distributed across \eqn{S}
#' remote servers. The average treatment effect (ATE) can be written as
#'
#' \deqn{
#' \widehat{ATE}\left(X; \widehat{\beta}_{ATE}\right)
#' =
#' \sum_{s=1}^{S}
#' \sum_{h=1}^{H}
#' \sum_{i=1}^{n_{s,h,d}}
#' w^{*}_{s,h,d,i}
#' \left(
#' \widehat{y}_{s,h,d=1,i}
#' -
#' \widehat{y}_{s,h,d=0,i}
#' \right)
#' }
#'
#' where \eqn{\widehat{y}_{ij} = Y_i \mid X_i, D_i = j} denotes the predicted outcome
#' for the \eqn{i}-th sample conditional on covariates \eqn{X_i} and treatment status
#' \eqn{j}. The corresponding weight is defined as
#'
#' \deqn{
#' w^{*}_{s,h,d,i}
#' =
#' \frac{n_{s,d}}{N_d}
#' \frac{n_{s,h,d}}{n_{s,d}}
#' \frac{w_{s,h,d,i}}{n_{s,h,d}}
#' }
#'
#' and depends on the number of samples within each stratum \eqn{h} and server
#' \eqn{s}, as well as on the matching procedure used to generate
#' \eqn{w_{s,h,d,i}}.
#'
#' The two innermost summations can be aggregated locally on each server prior to
#' communication with the central client. Consequently, only a weighted ATE from
#' each of the \eqn{S} remote servers is required to compute the global ATE estimate.
#'
#' @author Roy Gusinow
#'
#' @export
ds.avg_compute <- function(fit,
                           data,
                           treatment,
                           estimand = "ATE",
                           comparison = "difference",
                           error_type = NULL,
                           clusters = NULL,
                           eps = 1e-7,
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

  if (!is.character(treatment) || length(treatment) != 1) {
    stop("The 'treatment' argument must be a character string indicating a variable in the data.")
  }

  if (!estimand %in% c("ATE", "ATT", "ATC")) {
    stop("The 'estimand' must be one of: 'ATE', 'ATT', or 'ATC'.")
  }

  if (!comparison %in% c("difference", "lnriskratio", "lnoddsratio")) {
    stop("The 'comparison' must be one of: 'difference', 'lnriskratio' or 'lnoddsratio'.")
  }

  if (!error_type %in% c("HC0", "HC1", "HC2", "HC3", "const")) {
    stop("The 'error_type' must be one of 'HC0', 'HC1', 'HC2', 'HC3', or 'const'.")
  }

  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) {
    stop("The 'eps' argument must be a single positive numeric value.")
  }

  # fit$formula <- strip_weights_from_formula(fit$formula)
  avg_list <- list()

  # estimate
  estimate <- avg_compute_estimate(fit,
                                   data,
                                   treatment,
                                   estimand,
                                   comparison,
                                   datasources)
  avg_list$estimate <- estimate

  # standard error
  if (!is.null(error_type)){

    J <- ds.jacobian(fit,
                     data,
                     treatment,
                     estimand,
                     comparison,
                     eps,
                     datasources)

    vc <- ds.vcovCL(fit,
                    data,
                    error_type,
                    clusters,
                    datasources)
    se <- sqrt(t(J) %*% vc %*% J)

    avg_list$se <- se[1]
  }

  return(avg_list)
}

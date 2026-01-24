#' This function computes a heteroskedasticity-consistent (HC) sandwich variance-covariance matrix
#' for a GLM model fitted in a federated DataSHIELD setting. The sandwich estimator supports
#' HC0–HC3 types and optional cluster-robust inference.
#'
#' @param fit A list. Output from `ds.glm`, containing at least `VarCovMatrix`, `nsubs`, and `coefficients`.
#' @param data A character string. The name of the data frame on the server-side.
#' @param error_type A character string. Type of heteroskedasticity-consistent estimator to use. One of `"HC0"`, `"HC1"`, `"HC2"`, or `"HC3"`. Default is `"HC0"`.
#' @param clusters A numeric vector. Subclass labels (or other cluster IDs) used to aggregate estimating functions.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function uses `datashield.connections_find()` to detect active connections.
#'
#' @return A numeric matrix representing the robust variance-covariance matrix of the estimated model parameters.
#'
#' @details
#' The function constructs the sandwich estimator as \( B × M × B / n \), where:
#' \itemize{
#'   \item \( B \) is the scaled variance-covariance matrix ("bread") adjusted for dispersion.
#'   \item \( M \) is the cluster-aggregated estimating function matrix ("meat"), computed using `ds.estfun`.
#' }
#' This implementation supports cluster-robust corrections using the `meatCL_ds()` helper, and leverages
#' server-side assignments for computing hat matrix diagonals (`ds.hat_diag`) and score functions.
#'
#' The cluster-robust variance--covariance matrix is computed using the sandwich
#' estimator
#' \deqn{
#' \widehat{V} = \widehat{A}\,\widehat{B}\,\widehat{A}
#' }
#' where the ``bread'' estimator \eqn{\widehat{A}} is based on the empirical inverse
#' Hessian of the objective function \eqn{\Psi(y, x, \beta)}. The ``meat'' estimator
#' is given by
#'
#' \deqn{
#' \widehat{B}
#' =
#' \frac{1}{N}
#' \sum_{s=1}^{S}
#' \sum_{h=1}^{H}
#' \left(
#' \sum_{i=1}^{n_{s,h}}
#' \psi\!\left(
#' y_{s,h,i},
#' x_{s,h,i},
#' \widehat{\beta}_{ATE}
#' \right)
#' \right)
#' \left(
#' \sum_{i=1}^{n_{s,h}}
#' \psi\!\left(
#' y_{s,h,i},
#' x_{s,h,i},
#' \widehat{\beta}_{ATE}
#' \right)
#' \right)^{T}
#' }
#'
#' where \eqn{\psi} denotes the estimating function. For generalized linear models
#' implemented within DataSHIELD, \eqn{\psi} is equivalent to the evaluated score
#' function, i.e.\ the first derivative of the log-likelihood \eqn{\Psi}. In this
#' setting, the quantities communicated to the central client are products of
#' residuals and the design matrix, which are not recoverable at the individual
#' level once transmitted.
#'
#' When the number of clusters \eqn{H} is small, clusters residing on the same
#' server may be aggregated locally prior to communication with the central client.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item None, but uses `ds.dispersion`, `ds.hat_diag`, and `ds.estfun` to compute necessary components.
#' }
#'
#' @author Roy Gusinow
#' @export
ds.vcovCL <- function(fit,
                      data,
                      error_type = "HC0",
                      clusters,
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

  if (!error_type %in% c("HC0", "HC1", "HC2", "HC3", "const")) {
    stop("The 'error_type' must be one of 'HC0', 'HC1', 'HC2', 'HC3', or 'const'.")
  }

  # bread matrix
  dispersion <- ds.dispersion(fit,
                              data,
                              datasources)
  B <- fit$VarCovMatrix * fit$nsubs * dispersion

  # hat matrix diagonals assignment
  h_diag <- ds.hat_diag(fit,
                        data,
                        newobj_name = "h_diag",
                        datasources)
  # meat matrix
  if (!is.null(clusters)) {
    # when conducting cl errors, need to pass clusters to estfun
    # estfun returns rval unscaled list already aggregated by clusters
    rval_unscaled_list <- ds.estfun(fit,
                           data,
                           error_type,
                           h_diag = "h_diag",
                           using_clusters = TRUE,
                           datasources)

    score_eval <- do.call("rbind", rval_unscaled_list)

    M <- meatCL_ds(score_eval,
                   clusters,
                   error_type)
  }else{

    # when no cl errors, just get list of meats
    n <- fit$Nvalid # total number of samples
    k <- ncol(fit$VarCovMatrix) # number of parameters

    M_list <- ds.estfun(fit,
                           data,
                           error_type,
                           h_diag = "h_diag",
                           using_clusters = FALSE,
                           datasources)

    M <- Reduce(`+`, M_list)

    M <- meatHC_ds(M, error_type, n, k)
  }
  # rval_unscaled_list <- ds.estfun(fit,
  #                        data,
  #                        error_type,
  #                        h_diag = "h_diag",
  #                        datasources)
  #
  # score_eval <- do.call("rbind", rval_unscaled_list)

  # if (!is.null(clusters)) {
  #   M <- meatCL_ds(score_eval,
  #                  clusters,
  #                  error_type)
  # }else{
  #   M <- meatHC_ds(score_eval, error_type)
  # }

  # sand <- B %*% M %*% B / NROW(score_eval)
  sand <- B %*% M %*% B / fit$Nvalid

  return(sand)
}



#' This function performs propensity score matching in a federated DataSHIELD setting.
#' It estimates propensity scores using a logistic GLM, applies privacy-preserving noise,
#' and matches units based on the noisy scores using the `MatchIt` package logic on the client side.
#'
#' @param formula A formula object. Specifies the GLM used to estimate the propensity score. The treatment variable must appear on the left-hand side.
#' @param data A character string. The name of the data frame on the server-side.
#' @param sd A numeric value. The standard deviation of the Gaussian noise when using global noise
#' @param k An integer. The number of neighbors when using knn
#' @param newobj A character string. The name of the new object to create on each server containing the matched dataset. Default is `"matched_df"`.
#' @param datasources A list of `DSConnection-class` objects. If `NULL`, the function uses `datashield.connections_find()` to detect active connections.
#' @param ... Additional arguments passed to the underlying `MatchIt::matchit()` function on the client side (e.g., `method`, `distance`, `replace`, `m.order`).
#'
#' @return A list containing:
#' \describe{
#'   \item{subclass}{A vector of subclass labels for the matched samples.}
#' }
#'
#' @details
#' This function fits a logistic regression model using federated GLM (`ds.glm`) to estimate
#' the propensity score. It then adds noise according to the selected privacy scheme and
#' performs client-side matching using `MatchIt`. The final matched data, including weights,
#' is re-assigned to the servers via `matchitDS2`.
#'
#' **Server-side functions called**:
#' \itemize{
#'   \item `assign_IDDS`
#'   \item `matchitDS`
#'   \item `matchitDS2`
#' }
#'
#' @author Roy Gusinow
#'
#' @export
ds.matchit <- function(formula = NULL,
                       data = NULL,
                       sd = NULL,
                       k = NULL,
                       newobj = "matched_df",
                       datasources = NULL,
                       ...
){

  # look for DS connections
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }

  # if the argument 'data' is set, check that the data frame is defined (i.e. exists) on the server site
  if((is.null(data))){
    stop("Please provide a valid dataframe on the serverside")
    # defined <- isDefined(datasources, data)
  }

  # Ensure fit$formula is a formula; convert if not
  if (!inherits(formula, "formula")) {
    formula <- stats::as.formula(formula)
  }

  id_name <- "ID" # temp

  # generate unique ID to df
  cally <- call("assign_IDDS", data, id_name = id_name)
  temp <- DSI::datashield.assign.expr(datasources, data, cally)

  # generate propensity on the server
  fed_glm <- dsBaseClient::ds.glm(formula = formula,
                    data = data,
                    family = "binomial",
                    datasources = datasources)

  # gen prop scores on the server side - using fed coeffs
  coef_vec <- fed_glm$coefficients[, "Estimate"]
  names(coef_vec)[names(coef_vec) == "(Intercept)"] <- "Intercept"

  ds.genProp(formula = formula, coefficients = coef_vec, data = data, newobj = "distance", link = "logit", datasources = datasources)

  # delete weights col
  ds.delete_col(
    data = data,
    col_name = "distance",
    datasources = datasources
  )

  dsBaseClient::ds.cbind(x = c(data, "distance"),
           newobj = data,
           datasources = datasources)

  # begin matching
  # first call
  cally <- call("matchitDS", formula, data = data, distance = "distance", id_name = id_name, k, sd)
  noise_distance <- DSI::datashield.aggregate(datasources, cally)
  if (inherits(noise_distance, "character")){
    return(noise_distance)
  }

  # merge server prop frames and sort by ID
  prop_pool <- pool_servers(noise_distance)
  for (coln in all.vars(formula)){
    if(!(coln %in% colnames(prop_pool)))
    {
      prop_pool[, coln] <- 0
    }
  }

  # assign matched samples to data on servers
  # base_args <- list(
  #   formula  = formula,
  #   data     = subset(prop_pool, select = -c(distance)),
  #   distance = prop_pool$distance
  # )
  base_args <- list(
    formula  = formula,
    data     = prop_pool[ , setdiff(names(prop_pool), "distance"), drop = FALSE],
    distance = prop_pool[["distance"]]
  )
  dot_args <- list(...)
  if (dot_args$method == "optimal"){
    dot_args$m.order <- NULL
    dot_args$replace <- NULL
  }

  pooled.match.log <- do.call(MatchIt::matchit, c(base_args, dot_args))
  pooled.match <- MatchIt::match.data(pooled.match.log)

  # second call
  cally2 <- call("matchitDS2",
                 data,
                 pooled.match[, id_name],
                 pooled.match[, "distance"],
                 pooled.match[, "weights"],
                 id_name)
  result <- DSI::datashield.assign.expr(datasources, newobj, cally2)

  return(list(prop_pool = prop_pool,
              pooled.match = pooled.match,
              subclass = pooled.match$subclass))

  # return(list(subclass = pooled.match$subclass))
}

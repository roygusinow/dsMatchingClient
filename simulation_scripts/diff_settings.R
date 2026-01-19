# ------------------------------------------------------------------------------
# this script contains functions to compare different matching settings between
# federated and non-federated implementations (Figure 3c). Please note that the default privacy budget of 10^5 implies the central and DS implementations will slightly differ

# devtools::install("dsMatching")
# devtools::install("dsMatchingClient")

# Load all the required packages
library(dsBase)
library(dsBaseClient)
library(resourcer)
library(DSLite)
library(dsMatching)
library(dsMatchingClient)

library(ggplot2)
library(geometry)
library(dplyr)

library("MatchIt")
library(sandwich)
library("marginaleffects")


split_it <- function(df){
  df1 <- df[1:50,]
  df2 <- df[51:200,]
  df3 <- df[201:dim(df)[1],]

  return(list("df1" = df1, "df2" = df2, "df3" = df3))
}
compare_lists <- function(a, b) {
  # unique elements (if you care about duplicates, drop the unique() calls)
  a_u <- unique(a)
  b_u <- unique(b)

  # 1) intersection
  common   <- intersect(a_u, b_u)
  n_common <- length(common)

  # 2) “only in a” and “only in b”
  only_a   <- setdiff(a_u, b_u)
  only_b   <- setdiff(b_u, a_u)

  # 3) report
  cat("Number of shared elements:", n_common, "\n")
  if (n_common > 0) {
    cat("Shared elements:\n");   print(common)
  }

  if (length(only_a) > 0) {
    cat("In first list but not second:\n"); print(only_a)
  } else {
    cat("No elements unique to the first list.\n")
  }

  if (length(only_b) > 0) {
    cat("In second list but not first:\n"); print(only_b)
  } else {
    cat("No elements unique to the second list.\n")
  }
}
# ------------------------------------------------------------------------------
# non fed
# manually error compute

# ATE
estfun_man <- function(x,
                       type){
  # get evaluation of score vector (diff log likelihood) of

  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if(any(alias <- is.na(coef(x)))) xmat <- xmat[, !alias, drop = FALSE]
  wres <- as.vector(residuals(x, "working")) * weights(x, "working")
  dispersion <- if(substr(x$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  else sum(wres^2, na.rm = TRUE)/sum(weights(x, "working"), na.rm = TRUE)
  rval <- wres * xmat / dispersion
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL

  # temp addition
  # HC2/3
  if (type == "HC2" | type == "HC3"){
    # X <-  xmat
    h_diag <- hatvalues(x)

    if (type == "HC2"){
      rval <- rval / sqrt(1 - h_diag)
    }else if (type == "HC3"){
      rval <- rval / (1 - h_diag)
    }

  }

  return(rval)

}

bread_man <- function(model){

  if(!is.null(model$na.action)) class(model$na.action) <- "omit"
  sx <- summary(model)
  wres <- as.vector(residuals(model, "working")) * weights(model, "working")
  dispersion <- if(substr(model$family$family, 1L, 17L) %in% c("poisson", "binomial", "Negative Binomial")) 1
  else sum(wres^2)/sum(weights(model, "working"))

  print("whatever")
  print(as.vector(sum(sx$df[1L:2L])))

  return(sx$cov.unscaled * as.vector(sum(sx$df[1L:2L])) * dispersion)
}

meatCL_man <- function(model, clusters, type = "HC0"){


  score_eval <- estfun_man(model, type)
  # calculate the clustered with evals and the corresponding clusters
  n <- NROW(score_eval)
  k <- NCOL(score_eval)
  g <- length(unique(clusters)) # n of clusters

  # HC2/3
  if (type == "HC2" | type == "HC3"){
    score_eval <- sqrt((g - 1) / g) * score_eval
  }
  ## aggregate within cluster levels
  score_eval_agg <- apply(score_eval, 2L, rowsum, clusters)

  # outer product
  # adj <- 1
  adj <- g / (g - 1) # 1.002273
  rval <- adj * (t(score_eval_agg) %*% score_eval_agg)/n

  # HC0/1
  if (type == "HC0"){
    rval <- rval
  }else if (type == "HC1"){
    rval <- (n - 1)/(n - k) * rval
  }

  return(rval)

}

vcovCL_man <- function(model, clusters, type = "HC1"){
  # manually get vcov clustered -- currently constructed HC1 adj

  B <- bread_man(model)
  M <- meatCL_man(model, clusters, type)

  sand <- (1 / NROW(estfun_man(model, type))) * (B %*% M %*% B )

  return(sand)
}

get_avg_diff <- function(model,
                         treatment,
                         newdata){

  df_off <- newdata; df_on <- newdata;
  df_off[, treatment] <- 0; df_on[, treatment] <- 1;

  pred_treatment_off <- predict(model, df_off, type = "response")
  pred_treatment_on <- predict(model, df_on, type = "response")

  avg_diff <- mean(pred_treatment_on) - mean(pred_treatment_off)

  return(
    avg_diff
  )
}

jacobian_man <- function(model, treatment, newdata, eps = 1e-7){

  clean_avg <-  get_avg_diff(model,
                             treatment,
                             newdata)


  J <- c()
  for (i in 1:length(model$coefficients)){
    perturbed_model <- model

    fit_coef_perturbed <- perturbed_model$coefficients
    fit_coef_perturbed[i] <- fit_coef_perturbed[i] + eps
    perturbed_model$coefficients <- fit_coef_perturbed

    perturbed_avg <-  get_avg_diff(perturbed_model,
                                   treatment,
                                   newdata)


    J <- c(J, (clean_avg - perturbed_avg) / eps)
  }
  J <- as.matrix(J)

  return(J)
}

se_estimate <- function(model,
                        treatment,
                        clusters,
                        newdata,

                        type = "HC1",
                        eps = 1e-7){

  J <- jacobian_man(model,
                    treatment,
                    newdata = newdata,
                    eps)

  vc <- vcovCL_man(model,
                   clusters,
                   type)

  output <- sqrt(t(J) %*% vc %*% J)

  return(output)

}

avg_comparison_man <- function(model,
                               treatment,
                               clusters,
                               newdata,

                               type = "HC1",
                               eps = 1e-7){

  estimate <- get_avg_diff(model,
                           treatment,
                           newdata)
  se <- se_estimate(model,
                    treatment,
                    clusters,
                    newdata,

                    type,
                    eps)

  return(list(estimate = estimate,
              se = se))
}
# ------------------------------------------------------------------------------
# avg let's gooooo federated
get_avg_diff_sums_FL <- function(model,
                                 treatment,
                                 newdata){

  df_off <- newdata; df_on <- newdata;
  df_off[, treatment] <- 0; df_on[, treatment] <- 1;

  pred_treatment_off <- predict(model, df_off, type = "response")
  pred_treatment_on <- predict(model, df_on, type = "response")

  output_list <- c(sum_treated = sum(pred_treatment_on),
                   n_treat = length(pred_treatment_on),
                   sum_control = sum(pred_treatment_off),
                   n_control = length(pred_treatment_off))
  return(output_list)
}

get_avg_diff_means_FL <- function(model_list,
                                  treatment_list,
                                  newdata_list){


  sum_list <- mapply(get_avg_diff_sums_FL,
                     model = model_list,
                     treatment = treatment_list,
                     newdata = newdata_list,

                     SIMPLIFY = TRUE,
                     USE.NAMES = TRUE)

  summed_tab <- colSums(t(sum_list))
  ATE <- summed_tab["sum_treated"] / summed_tab["n_treat"] - summed_tab["sum_control"] / summed_tab["n_control"]

  return(ATE[[1]])
}

# error estimates federated
# jacobian
jacobian_man_FL <- function(model,
                            treatment_list,
                            newdata_list,
                            eps = 1e-7){

  model_list <- rep(list(model), length(newdata_list))
  clean_avg <- get_avg_diff_means_FL(model_list,
                                     treatment_list,
                                     newdata_list)

  J <- c()
  for (i in 1:length(model$coefficients)){
    perturbed_model <- model

    fit_coef_perturbed <- perturbed_model$coefficients
    fit_coef_perturbed[i] <- fit_coef_perturbed[i] + eps
    perturbed_model$coefficients <- fit_coef_perturbed
    print(perturbed_model$coefficients)

    perturbed_model_list <- rep(list(perturbed_model), length(newdata_list))
    perturbed_avg <- get_avg_diff_means_FL(perturbed_model_list,
                                           treatment_list,
                                           newdata_list)

    J <- c(J, (clean_avg - perturbed_avg) / eps)
  }
  J <- as.matrix(J)

  return(J)
}

# varcov
# bread
bread_FL <- function(
    global_res_df,
    global_cov_scaled){

  print("global_res_df")
  print(global_res_df)

  return(global_cov_scaled * global_res_df)
}
get_residuals <- function(y, pred){
  return(y - pred)
}
get_dispersion <-  function(model,
                            new_data,
                            family = "gaussian"){

  xmat <- model.matrix(formula(model$formula), new_data)

  outcome_var <- all.vars(formula(model$formula))[1] # get outcome var
  pred <- predict(model, new_data)
  wres <- get_residuals(new_data[, outcome_var], pred) * new_data[, "weights"] # weights from regressions model!

  dispersion <- if(family %in% c("poisson", "binomial", "Negative Binomial")) 1
  else sum(wres^2, na.rm = TRUE)

  return(dispersion)
}
get_summed_weight <- function(new_data){
  w <- new_data[, "weights"]
  return(sum(w))
}
get_H_diag <- function(model){
  out  <- hatvalues(model)
  return(out)
}
get_H_diag_man <- function(global_cov_mat,
                           xmat,
                           global_avg_vec,
                           global_n){

  # get hat diag using mahalanobis -- server side
  vec_diff <- xmat - rep(global_avg_vec, each = nrow(xmat))
  m_D2 <- vec_diff %*% solve(global_cov_mat) %*% t(vec_diff) # mahalanobis distance

  H_diag <- 1 / global_n + m_D2 / (global_n - 1)

  return(diag(H_diag))
}
get_H_list <- function(model,
                       newdata_list,
                       global_cov_mat,
                       global_avg_vec,
                       global_n){

  xmat_list <- list()
  for (i in 1:length(newdata_list)){
    a <- (model.matrix(formula(model$formula), newdata_list[[i]]))[, 2:dim(newdata_list)[2]]
    xmat_list[[i]] <- a
  }

  H_list <- mapply(get_H_diag_man,
                   global_cov_mat = global_cov_mat,
                   xmat = xmat_list,
                   global_avg_vec = global_avg_vec,
                   global_n = global_n,
                   SIMPLIFY = F)

  return(H_list)

}
# meat
get_rval_unscaled <- function(model,
                              new_data,
                              type,

                              global_cov_mat = NULL,
                              global_avg_vec = NULL,
                              global_n = NULL
){
  # get evaluation of score vector (diff log likelihood) of
  xmat <- model.matrix(formula(model$formula), new_data)

  outcome_var <- all.vars(formula(model$formula))[1] # get outcome var
  pred <- predict(model, new_data)
  wres <- get_residuals(new_data[, outcome_var], pred) * new_data[, "weights"] # weights from regressions model!
  rval_unscaled <- wres * xmat

  # HC2/3
  if (type == "HC2" | type == "HC3"){
    k <- dim(xmat)[2]
    h_diag <- get_H_diag_man(global_cov_mat,
                             xmat = xmat[, 2:k], # remove intercept
                             global_avg_vec,
                             global_n)


    if (type == "HC2"){
      rval_unscaled <- rval_unscaled / sqrt(1 - h_diag)
    }else if (type == "HC3"){
      rval_unscaled <- rval_unscaled / (1 - h_diag)
    }

  }

  return(rval_unscaled)

}
estman_FL <- function(model,
                      newdata_list,
                      type,

                      global_cov_mat = NULL,
                      global_avg_vec = NULL,
                      global_n = NULL){

  n_servers <- length(newdata_list)
  # get dispersion
  dispersion_list <- mapply(get_dispersion,
                            rep(list(model), n_servers),
                            newdata_list) # at server
  weight_summed_list <- lapply(newdata_list, get_summed_weight) # at server
  dispersion <- Reduce("+", dispersion_list) / Reduce("+", weight_summed_list)

  rval_unscaled_list <- mapply(get_rval_unscaled,
                               model = rep(list(model), n_servers),
                               new_data = newdata_list,
                               type = type,

                               global_cov_mat,
                               global_avg_vec,
                               global_n
  )
  rval_scaled_list <- mapply(function(x, y) x / y,
                             rval_unscaled_list,
                             rep(dispersion, n_servers))

  return(rval_scaled_list)
}
meatCL_ds <- function(score_eval,
                      clusters,
                      error_type){

  # calculate the clustered with evals and the corresponding clusters
  n <- NROW(score_eval)
  k <- NCOL(score_eval)
  g <- length(unique(clusters)) # n of clusters

  # HC2/3
  if (error_type == "HC2" | error_type == "HC3"){
    score_eval <- sqrt((g - 1) / g) * score_eval
  }
  ## aggregate within cluster levels
  score_eval_agg <- apply(score_eval, 2L, rowsum, clusters)

  # outer product
  # adj <- 1
  adj <- g / (g - 1)
  rval <- adj * (t(score_eval_agg) %*% score_eval_agg)/n

  # HC0/1
  if (error_type == "HC0"){
    rval <- rval
  }else if (error_type == "HC1"){
    rval <- (n - 1)/(n - k) * rval
  }

  return(rval)
}
meatCL_FL <- function(model,
                      newdata_list,
                      cluster_name,
                      type = "HC1",

                      global_cov_mat = NULL,
                      global_avg_vec = NULL,
                      global_n = NULL){

  # on the server side
  score_eval_list <- estman_FL(model,
                               newdata_list,
                               type,

                               global_cov_mat,
                               global_avg_vec,
                               global_n) # residuals * x at server

  score_eval <- do.call("rbind", score_eval_list)
  # calculate the clustered with evals and the corresponding clusters
  n <- NROW(score_eval)
  k <- NCOL(score_eval)

  clusters <- unlist(lapply(newdata_list, function(x) x[, cluster_name]))
  g <- length(unique(unlist(clusters))) # n of clusters

  # HC2/3
  if (type == "HC2" | type == "HC3"){
    score_eval <- sqrt((g - 1) / g) * score_eval
  }
  ## aggregate within cluster levels
  score_eval_agg <- apply(score_eval, 2L, rowsum, clusters)

  # outer product
  # adj <- (n - k)/(n - 1)
  adj <- g / (g - 1) # 1.002273
  rval <- adj * (t(score_eval_agg) %*% score_eval_agg)/n

  # HC0/1
  if (type == "HC0"){
    rval <- rval
  }else if (type == "HC1"){
    rval <- (n - 1)/(n - k) * rval
  }

  return(rval)
}
vcovCL_FL <- function(model,
                      newdata_list,
                      cluster_name,
                      type = "HC1"){
  # manually get vcov clustered -- currently constructed HC1 adj

  n_servers <- length(newdata_list)
  # get dispersion
  dispersion_list <- mapply(get_dispersion,
                            rep(list(model), n_servers),
                            newdata_list) # at server
  weight_summed_list <- lapply(newdata_list, get_summed_weight) # at server
  dispersion <- Reduce("+", dispersion_list) / Reduce("+", weight_summed_list)

  print("manual dispersion")
  print(dispersion)

  global_cov_scaled <- summary(model)$cov.unscaled * dispersion # unscaled * dispersion!
  global_res_df <- sum(summary(model)$df[1:2])

  B <- bread_FL(global_res_df,
                global_cov_scaled)

  # hat value params -- global vars
  xmat <- model.matrix(model)
  xmat <- xmat[, 2:dim(xmat)[2]]
  cov_mat <- cov(xmat)

  M <- meatCL_FL(model, newdata_list, cluster_name, type,

                 global_cov_mat = rep(list(cov_mat), n_servers),
                 global_avg_vec = rep(list(colMeans(xmat)), n_servers),
                 global_n = rep(nrow(xmat), n_servers))

  print(NROW(do.call("rbind", newdata_list)))
  sand <- B %*% M %*% B / NROW(do.call("rbind", newdata_list))

  return(sand)
}
se_estimate_FL <- function(model,
                           treatment_name,
                           cluster_name,
                           newdata_list,

                           type = "HC1",
                           eps = 1e-7){

  n_servers <- length(newdata_list)
  J <- jacobian_man_FL(model = model,
                       treatment_list = rep(treatment_name, n_servers),
                       newdata_list = newdata_list,
                       eps = eps)

  vc <- vcovCL_FL(model,
                  newdata_list,
                  cluster_name,
                  type)

  output <- sqrt(t(J) %*% vc %*% J)

  return(output)

}
avg_comparison_FL <- function(model,
                              treatment_name,
                              cluster_name,
                              newdata_list,

                              type = "HC1",
                              eps = 1e-7){

  n_servers <- length(newdata_list)
  estimate <- get_avg_diff_means_FL(model_list = rep(list(model), n_servers),
                                    treatment_list = rep(treatment_name, n_servers),
                                    newdata_list = newdata_list)
  se <- se_estimate_FL(model,
                       treatment_name,
                       cluster_name,
                       newdata_list,
                       type,
                       eps = eps)

  return(list(estimate = estimate,
              se = se))
}

# ------------------------------------------------------------------------------
# check if the results are similar
check_equivalence <- function(
    data_df,
    params_list,
    comparison_tol = 1e-7
){

  # get params
  matching_form <- params_list[["matching_form"]]
  est_form <- params_list[["est_form"]]
  family <- params_list[["family"]]
  meth <- params_list[["meth"]]
  m.order <- params_list[["m.order"]]
  ratio <- params_list[["ratio"]]
  replace <- params_list[["replace"]]
  tol <- params_list[["tol"]]
  discard <- params_list[["discard"]]
  eps <- params_list[["eps"]]
  correction_type <- params_list[["correction_type"]]

  # localised matching -- handling of bug
  if (meth == "optimal" && (is.null(m.order) || !is.null(replace))){
    norm_match <- matchit(formula = formula(matching_form),
                          data = data_df,

                          method = meth,
                          # m.order = m.order,
                          ratio = ratio,
                          # replace = replace,
                          tol = tol,
                          discard = discard)
  }else{
    norm_match <- matchit(formula = formula(matching_form),
                          data = data_df,

                          method = meth,
                          m.order = m.order,
                          ratio = ratio,
                          replace = replace,
                          tol = tol,
                          discard = discard)
  }
  md <- match.data(norm_match)

  # normal error compute
  fit_non_FL_standard <- glm(formula = est_form,
                             data = md,
                             family = family)
  # md$subclass <- factor(1:length(md$subclass))
  clusters <- md$subclass
  vc_auto <-vcovCL(fit_non_FL_standard, clusters, type = correction_type)
  avg_non_FL_standard <- avg_comparisons(fit_non_FL_standard , variables = "A",
                                         newdata = md,
                                         wts = "weights",
                                         vcov = vc_auto,

                                         eps = eps,
                                         numderiv = "fdforward")

  print("Standard Package Results")
  print(avg_non_FL_standard)


  # manual error computed (non federated)
  avg_non_FL_man <- avg_comparison_man(fit_non_FL_standard,
                                       "A",
                                       clusters,
                                       md,

                                       type = correction_type,
                                       eps = eps
  )
  print("Man Results")
  print(avg_non_FL_man)


  if (abs(avg_non_FL_man$se - avg_non_FL_standard$std.error) > comparison_tol){
    print("WARNING: Standard Package and Manual implementation error estimates are not equivalent")

    print("Checking vcov")

    print("Non FL vcov")
    print(vc_auto)

    print("Manual vcov")
    vc_man <- vcovCL_man(fit_non_FL_standard,
                         clusters,
                         correction_type)
    print(vc_man)

  }else{
    print("SUCCESS: Standard Package and Manual implementation error estimates are equivalent up to threshold")
  }

  # manual error computed (federated)
  sp <- split_it(md)
  newdata_list <- list(sp$df1,
                       sp$df2,
                       sp$df3)
  n_servers <- length(newdata_list)
  avg_FL_local_std <- avg_comparison_FL(fit_non_FL_standard,
                                        treatment_name = "A",
                                        cluster_name = "subclass",
                                        newdata_list,

                                        type = correction_type,
                                        eps = eps
  )
  if (abs(avg_non_FL_man$se - avg_FL_local_std$se) > comparison_tol){
    print("WARNING: Manual and FL are not equivalent")

    print("Checking vcov")

    print("Manual vcov")
    vc_man <- vcovCL_man(fit_non_FL_standard,
                         clusters,
                         correction_type)
    print(vc_man)


    print("FL vcov")
    vc_FL <- vcovCL_FL(fit_non_FL_standard,
                       newdata_list,
                       cluster_name = "subclass",
                       correction_type)
    print(vc_FL)

  }else{
    print("SUCCESS: Manual and FL are equivalent")
  }

  # datashield test
  fed_match <- ds.matchit(form = formula(matching_form),
                          data = "DST",
                          privacy = "norm",
                          sd = 0,
                          # k = k,
                          newobj = "matched_pooled",

                          method = meth,
                          m.order = m.order,
                          ratio = ratio,
                          replace = replace,
                          tol = tol,
                          discard = discard)

  ds_fit <- ds.glm(formula = formula(est_form),
                   data = "matched_pooled",
                   family = family,
                   viewVarCov = T)

  ds_ate <- ds.avg_compute(fit = ds_fit,
                           data = "matched_pooled",
                           treatment = "A",
                           # avg_type = "ATE",

                           error_type = correction_type,
                           clusters = fed_match$subclass,
                           eps = eps
  )
  if (abs(ds_ate$se - avg_non_FL_man$se) > comparison_tol){
    print("WARNING: Datashield and Manual are not equivalent")
    print(abs(ds_ate$se - avg_non_FL_man$se))
    out <- "Failure"

    matched_pooled <- getDSLiteData(connections, "matched_pooled")
    id_v <- c()
    for (i in 1:length(matched_pooled)) {
      id_v <- c(id_v, rownames(matched_pooled[[i]]))
    }
    compare_lists(rownames(md), id_v)

    print("Checking vcov")

    vc_FL <- vcovCL_FL(fit_non_FL_standard,
                       newdata_list,
                       cluster_name = "subclass",
                       correction_type)
    print("FL vcov")
    print(vc_FL)


    print("DataSHIELD vcov")
    vc_ds <- ds.vcovCL(ds_fit,
                       "matched_pooled",
                       error_type = correction_type,
                       clusters = fed_match$subclass)
    print(vc_ds)

    print("Checking subclasses")
    print("fed match subclasses")
    print(fed_match$subclass)
    print("md subclasses")
    print(md$subclass)
    print("matched_pooled  fed")
    print(matched_pooled$subclass)

    print("Checking bread values")
    print("bread (standard bread)")
    bread_man <- bread_man(fit_non_FL_standard)
    print(bread_man)

    print("datashield bread")
    dispersion <- ds.dispersion(ds_fit,
                                data = "matched_pooled",
                                datasources = connections)
    print("datashield dispersion")
    print(dispersion)

    # res_df <- ds_fit$nsubs - NROW(ds_fit$coefficients)
    bread_ds <- ds_fit$VarCovMatrix * ds_fit$nsubs * dispersion
    print(bread_ds)
    print("compare bread")
    print(max(abs(bread_ds - bread_man)))

    print("Checking meat values")
    print("est function")
    EST_MAN <- estfun_man(fit_non_FL_standard, correction_type)
    # print(EST_MAN)

    print("meat (standard meat)")
    meat_man <- meatCL_man(fit_non_FL_standard, md$subclass, type = correction_type)
    print(meat(fit_non_FL_standard))
    print(meat_man)

    h_diag <- ds.hat_diag(ds_fit,
                          data = "matched_pooled",
                          newobj_name = "h_diag",
                          datasources = connections)
    # meat matrix
    estfun_FL <- ds.estfun(ds_fit,
                           data = "matched_pooled",
                           error_type = correction_type,
                           h_diag = "h_diag",
                           datasources  = connections)
    score_eval <- do.call("rbind", estfun_FL)
    # print(score_eval)

    print("compare est value max dif")
    print(max(abs(EST_MAN - score_eval)))
    print("compare clusers")


    print("datashield meat")
    meat_ds <- meatCL_ds(score_eval,
                         fed_match$subclass,
                         # md$subclass,
                         error_type = correction_type)
    print(meat_ds)
    print("compare meat")
    print(max(abs(meat_ds - meat_man)))

    print("try use real bread for vcov?")
    vcov_real_bread <- bread_man %*% meat_ds %*% bread_man / NROW(score_eval)
    print(vcov_real_bread)
    print("max difs between vcov and real_bread_vcov")
    print(max(abs(vcov_real_bread - vc_FL)))

    print("sanity check -- let's do bread again")
    vcov_norm_bread <- bread_ds %*% meat_ds %*% bread_ds / NROW(score_eval)
    print(vcov_norm_bread)
    print("max difs between vcov and real_bread_vcov")
    print(max(abs(vcov_norm_bread - vc_FL)))

    print("compare normal covs")
    global_cov_mat <- summary(fit_non_FL_standard)$cov.unscaled
    print("normal cov unscaled")
    print(global_cov_mat)
    print("datashield cov unscaled")
    print(ds_fit$VarCovMatrix)

    print("print cov dif")
    print(max(abs(global_cov_mat - ds_fit$VarCovMatrix)))

    print("Checking jacobian")
    print("FL jac")
    n_servers <- length(newdata_list)
    J_man_FL <- jacobian_man_FL(model = fit_non_FL_standard,
                                treatment_list = rep("treat", n_servers),
                                newdata_list = newdata_list,
                                eps = eps
    )
    print(J_man_FL)

    print("DataSHIELD jac")
    J_ds <- ds.jacobian(ds_fit,
                        data = "matched_pooled",
                        treatment = "A",
                        estimand = "ATE",
                        eps = eps,
                        datasources = connections
    )
    print(J_ds)
    print("compare jacobian")
    max_diff_jac <- max(abs(J_ds - J_man_FL))
    print(max_diff_jac)

    print("use real bread vcov for complete error")
    se <- sqrt(t(J_ds) %*% vcov_real_bread %*% J_ds)
    print(se)

    print("sanity check use norm bread vcov for complete error")
    se <- sqrt(t(J_ds) %*% vcov_norm_bread %*% J_ds)
    print(se)
  }
  else{
    print("SUCCESS: Datashield and Manual are equivalent")
    out <- "Success"
  }

  print("Compare Estimates")
  print(avg_non_FL_standard$estimate)
  print(avg_non_FL_man$estimate)
  print(avg_FL_local_std$estimate)
  print(ds_ate$estimate)

  print("Compare Error Estimates")
  print(avg_non_FL_standard$std.error)
  print(avg_non_FL_man$se)
  print(avg_FL_local_std$se)
  print(ds_ate$se)

  out <- list(
    message = out,
    avg_non_FL_standard = avg_non_FL_standard,
    avg_non_FL_man = avg_non_FL_man,
    avg_FL_local_std = avg_FL_local_std,
    ds_ate = ds_ate
  )
  return(out)
}

# ------------------------------------------------------------------------------
# run it

# set seed
seed_num <- 741
set.seed(seed_num)

# real data
data_path <- "data/austin_simulated_data/sim_data.csv"

read_sim <- function(path) {
  df <- read.csv2(path, stringsAsFactors = FALSE)
  if (!("ID" %in% names(df))) df$ID <- seq_len(nrow(df))
  df
}
sim_df <- read_sim(data_path)

# remove duplicates
all_rel_vars <- " A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9"
vars <- all.vars(formula(all_rel_vars))
sim_df <- sim_df[!duplicated(sim_df[, vars]), ]

sp <- split_it(sim_df)
sim_df1 <- sp$df1
sim_df2 <- sp$df2
sim_df3 <- sp$df3

# # log into ds
options(datashield.errors.print = TRUE)
cfg <- DSLite::defaultDSConfiguration(include=c("dsBase", "dsMatching"))
cfg$Options$default.datashield.privacyControlLevel <- "\"permissive\""
cfg$Options$default.nfilter.glm <- "1"
cfg$Options$default.nfilter.noise <- "-1"

# Create virtualized server with DSLite, assign everything needed on it
dslite.server <- newDSLiteServer(
  tables=list(
    sim_df1=sim_df1, sim_df2=sim_df2, sim_df3=sim_df3),
  config = cfg
)

builder <- DSI::newDSLoginBuilder()
builder$append(server = "server1", url = "dslite.server", table = "sim_df1", driver = "DSLiteDriver")
builder$append(server = "server2", url = "dslite.server", table = "sim_df2", driver = "DSLiteDriver")
builder$append(server = "server3", url = "dslite.server", table = "sim_df3", driver = "DSLiteDriver")

logindata.dslite <- builder$build()

# Login to the virtualized server
connections <- datashield.login(logindata.dslite, assign=T)

tab_symb <- "DST"
DSI::datashield.assign.table(conns = connections, symbol = tab_symb, table = c("sim_df1","sim_df2", "sim_df3"))

dslite.server$aggregateMethod("print", function(x){ print(x) })

match_form  <- "A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9"
effect_form <- "Y ~ A + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9"

# list of param checks
ll_param_checks <- list(
  # nearest
  list(
    matching_form = match_form,
    est_form = effect_form,
    family = "gaussian",
    meth = "nearest",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE,
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC0"
  ),
  # ratio
  list(
    matching_form = match_form,
    est_form = effect_form,
    family = "gaussian",
    meth = "nearest",
    m.order = "smallest",
    ratio = 2,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC0"
  ),
  # discard
  list(
    matching_form = match_form,
    est_form = effect_form,
    family = "gaussian",
    meth = "nearest",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "none",

    eps = 1e-4, # delta method
    correction_type = "HC0"
  ),
  # hc1
  list(
    matching_form = match_form,
    est_form = effect_form,
    family = "gaussian",
    meth = "optimal",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC1"
  ),
  # hc2
  list(
    matching_form = match_form,
    est_form = effect_form,
    family = "gaussian",
    meth = "optimal",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC2"
  ),
  # hc3
  list(
    matching_form = match_form,
    est_form = effect_form,
    family = "gaussian",
    meth = "optimal",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC3"
  ),
  # optimal
  list(
    matching_form = match_form,
    est_form = effect_form,
    family = "gaussian",
    meth = "optimal",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC0"
  )

  # nearest - logistic/binomial
  # list(
  #   matching_form = "treat ~ age + educ + re74 + re75",
  #   # est_form = "re78 ~ treat + age + re75 + educ",
  #   # family = "gaussian",
  #
  #   est_form = "married ~ treat + age + re75 + educ",
  #   family = "binomial",
  #
  #   meth = "nearest",
  #   m.order = "smallest",
  #   ratio = 1,
  #   replace = FALSE,
  #   tol = 1e-7,
  #   discard = "both",
  #
  #   eps = 1e-4, # delta method
  #   correction_type = "HC0"
  # )

)
param_list_names = c("Nearest", "2:1 Ratio", "No Discarding", "HC1 Error", "HC2 Error", "HC3 Error", "Optimal Matching")

# comparison_tol <- 1
comparison_tol <- 10e-6

suc_list <- vector("list", length(ll_param_checks))

for (i in seq_along(ll_param_checks)) {
  message("Running test ", i)
  suc_list[[i]] <- check_equivalence(
    data_df = sim_df,
    params_list = ll_param_checks[[i]],
    comparison_tol = comparison_tol
  )
}

# save
diff_settings_df <- data.frame(
  type = character(),
  method = character(),
  ate_est = numeric(),
  ate_std_err = numeric()
)
for (i in 1:length(suc_list)){

  res <- suc_list[[i]]

  total_df <- data.frame(
    type = character(),
    method = character(),
    ate_est = numeric(),
    ate_std_err = numeric()
  )
  for (j in 2:length(res)) {
    temp_df <- data.frame(
      type = param_list_names[i],
      method = names(res)[j],
      ate_est = res[[j]]$estimate,
      ate_std_err = ifelse(j == 2, res[[j]]$std.error, res[[j]]$se)
    )

    total_df <- rbind(total_df, temp_df)
  }
  diff_settings_df <- rbind(diff_settings_df, total_df)

}

# write.csv2(diff_settings_df, "results/data_4_plotting/diff_settings.csv") # works!


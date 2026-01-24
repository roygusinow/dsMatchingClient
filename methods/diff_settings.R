# Script to recreate Figure 3c in paper. PLease note that results aren't
# exact as there is some randomness in matching (non-infinite privacy budget)

# Load all the required packages
library(dsBase)
library(dsBaseClient)
library(resourcer)
library(DSLite)
library(dsMatching)
library(dsMatchingClient)

library(dplyr)

library("MatchIt")
library(sandwich)
library("marginaleffects")

source("simulation.R")


split_it <- function(df){
  df1 <- df[1:50,]
  df2 <- df[51:200,]
  df3 <- df[201:dim(df)[1],]

  return(list("df1" = df1, "df2" = df2, "df3" = df3))
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
  comparison_type <- params_list[["comparison_type"]]

  if (comparison_type == "lnriskratio"){
    comparison_type_pack <- "lnratioavg"
  }else if (comparison_type == "lnoddsratio"){
    comparison_type_pack <- "lnoravg"
  }else{
    comparison_type_pack <- comparison_type
  }

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

  clusters <- md$subclass
  vc_auto <-vcovCL(fit_non_FL_standard, clusters, type = correction_type)
  avg_non_FL_standard <- avg_comparisons(fit_non_FL_standard , variables = "A",
                                         newdata = md,
                                         wts = "weights",
                                         vcov = vc_auto,
                                         comparison = comparison_type_pack,
                                         eps = eps,
                                         numderiv = "fdforward")

  # datashield test
  fed_match <- ds.matchit(form = formula(matching_form),
                          data = "DST",
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
                           estimand = "ATE",
                           comparison = comparison_type,
                           error_type = correction_type,
                           clusters = fed_match$subclass,
                           eps = eps
  )

  out <- list(
    avg_non_FL_standard = avg_non_FL_standard,
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
sim_df <- gen_data_correct(n = 1000, seed = seed_num)

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
    correction_type = "HC0",
    comparison_type = "difference"
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
    correction_type = "HC0",
    comparison_type = "difference"
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
    correction_type = "HC0",
    comparison_type = "difference"
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
    correction_type = "HC1",
    comparison_type = "difference"
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
    correction_type = "HC0",
    comparison_type = "difference"
  )

)
param_list_names = c("Nearest", "2:1 Ratio", "No Discarding", "HC1 Error", "Optimal Matching")

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
    type = rep(param_list_names[i], 2),
    method = c("Non-FL Standard", "DataSHIELD"),
    ate_est = c(res$avg_non_FL_standard$estimate, res$ds_ate$estimate),
    ate_std_err = c(res$avg_non_FL_standard$std.error, res$ds_ate$std.error)
  )
  diff_settings_df <- rbind(diff_settings_df, total_df)
}
print(diff_settings_df)


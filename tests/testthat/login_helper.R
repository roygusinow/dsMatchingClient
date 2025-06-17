create_test_data <- function(seed) {

  set.seed(seed)

  # real data
  data("lalonde", package = "MatchIt");

  # remove duplicates
  all_rel_vars <- "treat ~ age + educ + married + nodegree + re74 + re75"
  vars <- all.vars(formula(all_rel_vars))
  lalonde <- lalonde[!duplicated(lalonde[, vars]), ]

  shuffled_lalonde <- lalonde[sample(1:nrow(lalonde)), ]
  sp <- split_it(shuffled_lalonde)
  data_server_1 <- sp$df1
  data_server_2 <- sp$df2
  data_server_3 <- sp$df3

  # log into ds
  options(datashield.errors.print = TRUE)
  cfg <- DSLite::defaultDSConfiguration(include=c("dsBase", "dsMatching"))
  cfg$Options$default.datashield.privacyControlLevel <- "\"permissive\""
  cfg$Options$default.nfilter.glm <- "1"
  cfg$Options$default.nfilter.noise <- "-1"

  # Create virtualized server with DSLite, assign everything needed on it
  dslite.server <<- DSLite::newDSLiteServer(
    tables=list(
      data_server_1 = data_server_1,
      data_server_2 = data_server_2,
      data_server_3 = data_server_3),
    config = cfg
  )

  builder <- DSI::newDSLoginBuilder()
  builder$append(server = "server1", url = "dslite.server", table = "data_server_1", driver = "DSLiteDriver")
  builder$append(server = "server2", url = "dslite.server", table = "data_server_2", driver = "DSLiteDriver")
  builder$append(server = "server3", url = "dslite.server", table = "data_server_3", driver = "DSLiteDriver")

  logindata.dslite <- builder$build()

  # Login to the virtualized server
  connections <- DSI::datashield.login(logindata.dslite, assign=T)

  DSI::datashield.assign.table(conns = connections, symbol = "DST", table = c("data_server_1","data_server_2", "data_server_3"))

  return(list("connections" = connections, "data" = shuffled_lalonde))
}

create_test_data_one_server <- function(seed) {

  set.seed(seed)

  # real data
  data("lalonde", package = "MatchIt");

  # remove duplicates
  all_rel_vars <- "treat ~ age + educ + married + nodegree + re74 + re75"
  vars <- all.vars(formula(all_rel_vars))
  lalonde <- lalonde[!duplicated(lalonde[, vars]), ]

  shuffled_lalonde <- lalonde[sample(1:nrow(lalonde)), ]
  sp <- split_it(shuffled_lalonde)
  data_server_1 <- sp$df1
  data_server_2 <- sp$df2
  data_server_3 <- sp$df3

  # log into ds
  options(datashield.errors.print = TRUE)
  cfg <- DSLite::defaultDSConfiguration(include=c("dsBase", "dsMatching"))
  cfg$Options$default.datashield.privacyControlLevel <- "\"permissive\""
  cfg$Options$default.nfilter.glm <- "1"

  # Create virtualized server with DSLite, assign everything needed on it
  dslite.server <<- DSLite::newDSLiteServer(
    tables=list(
      data_server_1 = data_server_1),
    config = cfg
  )

  builder <- DSI::newDSLoginBuilder()
  builder$append(server = "server1", url = "dslite.server", table = "data_server_1", driver = "DSLiteDriver")

  logindata.dslite <- builder$build()

  # Login to the virtualized server
  connections <- DSI::datashield.login(logindata.dslite, assign=T)

  DSI::datashield.assign.table(conns = connections, symbol = "DST", table = c("data_server_1","data_server_2", "data_server_3"))

  return(list("connections" = connections, "data" = shuffled_lalonde))
}

split_it <- function(df){
  df1 <- df[1:50,]
  df2 <- df[51:200,]
  df3 <- df[201:dim(df)[1],]

  return(list("df1" = df1, "df2" = df2, "df3" = df3))
}


local_matching <- function(
    data_df,
    params_list
){

  # get params
  matching_form <- params_list[["matching_form"]]
  est_form <- params_list[["est_form"]]
  meth <- params_list[["meth"]]
  m.order <- params_list[["m.order"]]
  ratio <- params_list[["ratio"]]
  replace <- params_list[["replace"]]
  tol <- params_list[["tol"]]
  discard <- params_list[["discard"]]
  eps <- params_list[["eps"]]
  correction_type <- params_list[["correction_type"]]
  n_subclass <- params_list[["n_subclass"]]
  estimand <- params_list[["estimand"]]


  # localised matching
  norm_match <- MatchIt::matchit(formula = matching_form,
                                 data = data_df,
                                 subclass = n_subclass,

                                 method = meth,
                                 m.order = m.order,
                                 ratio = ratio,
                                 replace = replace,
                                 tol = tol,
                                 discard = discard)

  return(norm_match)
}

avg_local_matching <- function(
    data_df,
    params_list
){

  # get params
  matching_form <- params_list[["matching_form"]]
  est_form <- params_list[["est_form"]]
  meth <- params_list[["meth"]]
  m.order <- params_list[["m.order"]]
  ratio <- params_list[["ratio"]]
  replace <- params_list[["replace"]]
  tol <- params_list[["tol"]]
  discard <- params_list[["discard"]]
  eps <- params_list[["eps"]]
  correction_type <- params_list[["correction_type"]]
  n_subclass <- params_list[["n_subclass"]]
  estimand <- params_list[["estimand"]]
  treatment <- params_list[["treatment"]]

  # localised matching
  norm_match <- local_matching(data_df, params_list)

  md <- MatchIt::match.data(norm_match)

  # normal error compute
  fit_non_FL_standard <- stats::glm(formula = est_form,
                             data = md,
                             family = gaussian)
  vc_auto <- sandwich::vcovCL(fit_non_FL_standard, md$subclass, type = correction_type)

  if (estimand == "ATE") {
    dt <- md
  } else if (estimand == "ATT") {
    dt <- md[md$treat == 1, ]
  } else if (estimand == "ATC") {
    dt <- md[md$treat == 0, ]
  } else {
    stop("Unknown estimand type")
  }
  avg_non_FL_standard <- marginaleffects::avg_comparisons(fit_non_FL_standard,
                                         variables = treatment,
                                         newdata = dt,
                                         wts = "weights",
                                         vcov = vc_auto,

                                         eps = eps,
                                         numderiv = "fdforward")
  return(avg_non_FL_standard)
}

local_matching_subclass <- function(
    data_df,
    params_list
){

  # get params
  matching_form <- params_list[["matching_form"]]
  est_form <- params_list[["est_form"]]
  meth <- params_list[["meth"]]
  m.order <- params_list[["m.order"]]
  ratio <- params_list[["ratio"]]
  replace <- params_list[["replace"]]
  tol <- params_list[["tol"]]
  discard <- params_list[["discard"]]
  eps <- params_list[["eps"]]
  correction_type <- params_list[["correction_type"]]
  n_subclass <- params_list[["n_subclass"]]
  estimand <- params_list[["estimand"]]

  norm_match <- MatchIt::matchit(formula = matching_form,
                                 data = data_df,
                                 subclass = n_subclass,
                                 estimand = estimand,

                                 method = meth,
                                 m.order = m.order,
                                 ratio = ratio,
                                 replace = replace,
                                 tol = tol,
                                 discard = discard)

  md <- MatchIt::match.data(norm_match)

  return(norm_match)
}

ds_matching <- function(
    connections,
    params_list
){

  # get params
  matching_form <- params_list[["matching_form"]]
  est_form <- params_list[["est_form"]]
  meth <- params_list[["meth"]]
  m.order <- params_list[["m.order"]]
  ratio <- params_list[["ratio"]]
  replace <- params_list[["replace"]]
  tol <- params_list[["tol"]]
  discard <- params_list[["discard"]]
  eps <- params_list[["eps"]]
  correction_type <- params_list[["correction_type"]]
  estimand <- params_list[["estimand"]]
  treatment <- params_list[["treatment"]]

  # datashield test
  fed_match <- ds.matchit(formula = matching_form,
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
                          discard = discard,
                          datasources = connections)

  ds_fit <- dsBaseClient::ds.glm(formula = est_form,
                   data = "matched_pooled",
                   family = "gaussian",
                   viewVarCov = T,
                   datasources = connections)

  ds_ate <- ds.avg_compute(fit = ds_fit,
                           data = "matched_pooled",
                           treatment = treatment,
                           estimand = estimand,

                           error_type = correction_type,
                           clusters = fed_match$subclass,
                           eps = eps,
                           datasources = connections
  )
  return(ds_ate)
}

ds_matching_subclass <- function(
    connections,
    params_list
){

  # get params
  matching_form <- params_list[["matching_form"]]
  est_form <- params_list[["est_form"]]
  meth <- params_list[["meth"]]
  m.order <- params_list[["m.order"]]
  ratio <- params_list[["ratio"]]
  replace <- params_list[["replace"]]
  tol <- params_list[["tol"]]
  discard <- params_list[["discard"]]
  eps <- params_list[["eps"]]
  correction_type <- params_list[["correction_type"]]
  estimand <- params_list[["estimand"]]
  n_subclass <- params_list[["n_subclass"]]

  # datashield test
  fed_match <- ds.matchit_subclass(formula = matching_form,
                                   data = "DST",
                                   subclass = n_subclass,
                                   estimand = estimand, # working
                                   len = 600,
                                   newobj = "matched_pooled",
                                   datasources = connections)

  # ds_fit <- ds.glm(formula = formula(est_form),
  #                  data = "matched_pooled",
  #                  family = "gaussian",
  #                  viewVarCov = T,
  #                  weights = "weights",
  #                  datasources = connections)

  # not clear how to confirm correctnesss..
  # ds_ate <- ds.avg_compute(fit = ds_fit,
  #                          data = "matched_pooled",
  #                          treatment = "treat",
  #                          estimand = estimand,
  #                          clusters = fed_match$subclass,
  #                          error_type = correction_type,
  #                          eps = eps,
  #                          datasources = connections
  # )
  return(fed_match)
}



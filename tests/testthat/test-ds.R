# load file
source("login_helper.R")

tolerance <- 1e-3

# test that
test_that("Basic DataSHIELD Check", {
  data_object <- create_test_data(seed = 12345)
  connections <- data_object$connections
  test_data <- data_object$data

  mean_age <- dsBaseClient::ds.mean(x = "DST$age", type = "combine", datasources = connections)
  mean_age <- as.numeric(mean_age$Global.Mean[1])

  expect_equal(mean(test_data$age), mean_age, )
})

test_that(
  "-- Basic Summary test (Compare Means Treated) --", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <- list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = re78 ~ treat + age + re75 + educ,
      meth = "nearest",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE,
      tol = 1e-7,
      discard = "both",

      eps = 1e-4, # delta method
      correction_type = "HC0",
      estimand = "ATE",
      treatment = "treat"
    )

    # normal matching
    norm_match <- MatchIt::matchit(formula = params_list[["matching_form"]],
                          data = test_data,

                          method = params_list[["meth"]],
                          m.order = params_list[["m.order"]],
                          ratio = params_list[["ratio"]],
                          replace = params_list[["replace"]],
                          tol = params_list[["tol"]],
                          discard = params_list[["discard"]])
    # summary of unmatched/matched results
    norm_summary <- summary(norm_match, un = T, improvement = T)

    # federated matching
    fed_match <- ds.matchit(form = params_list[["matching_form"]],
                            data = "DST",
                            privacy = "norm",
                            sd = 0,
                            newobj = "matched_pooled",

                            method = params_list[["meth"]],
                            m.order = params_list[["m.order"]],
                            ratio = params_list[["ratio"]],
                            replace = params_list[["replace"]],
                            tol = params_list[["tol"]],
                            discard = params_list[["discard"]],
                            datasources = connections)

    # federated summary
    fed_summary.combined <- ds.match_summary(unmatched_obj = "DST",
                                             matched_obj = "matched_pooled",
                                             type = "combined",
                                             bin_num = 429,
                                             treatment = params_list[["treatment"]],
                                             datasources = connections)

    vars_comparison <- all.vars(params_list[["matching_form"]])
    vars_comparison <- vars_comparison[vars_comparison != params_list[["treatment"]]]
    vars_comparison <- c("distance", vars_comparison)

    # select rows with same row names
    summed_comparison <- fed_summary.combined$`Summary of Balance for Matched Data:`[vars_comparison, ]
    compared_vec <- summed_comparison$Means.Treated - norm_summary$sum.matched[, "Means Treated"]

    expect_equal(0, max(compared_vec), tolerance=tolerance)
  })

test_that(
  "-- Love Plotting test --", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <- list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = re78 ~ treat + age + re75 + educ,
      meth = "nearest",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE,
      tol = 1e-7,
      discard = "both",

      eps = 1e-4, # delta method
      correction_type = "HC0",
      estimand = "ATE",
      treatment = "treat"
    )

    fed_match <- ds.matchit(form = params_list[["matching_form"]],
                            data = "DST",
                            privacy = "norm",
                            sd = 0,
                            newobj = "matched_pooled",

                            method = params_list[["meth"]],
                            m.order = params_list[["m.order"]],
                            ratio = params_list[["ratio"]],
                            replace = params_list[["replace"]],
                            tol = params_list[["tol"]],
                            discard = params_list[["discard"]],
                            datasources = connections)

    # federated summary
    fed_summary.combined <- ds.match_summary(unmatched_obj = "DST",
                                             matched_obj = "matched_pooled",
                                             type = "combined",
                                             bin_num = 429,
                                             treatment = params_list[["treatment"]],
                                             datasources = connections)
    p <- love.plot(fed_summary.combined)

    expect_s3_class(p, "gg")
  })

test_that(
  "-- Jitter Plotting test --", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <- list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = re78 ~ treat + age + re75 + educ,
      meth = "nearest",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE,
      tol = 1e-7,
      discard = "both",

      eps = 1e-4, # delta method
      correction_type = "HC0",
      estimand = "ATE",
      treatment = "treat"
    )

    fed_match <- ds.matchit(form = params_list[["matching_form"]],
                            data = "DST",
                            privacy = "norm",
                            sd = 0,
                            newobj = "matched_pooled",

                            method = params_list[["meth"]],
                            m.order = params_list[["m.order"]],
                            ratio = params_list[["ratio"]],
                            replace = params_list[["replace"]],
                            tol = params_list[["tol"]],
                            discard = params_list[["discard"]],
                            datasources = connections)

    p <- jitter.plot(fed_match)

    expect_s3_class(p, "gtable")
  })

test_that(
  "-- QQ Plotting test --", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <- list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = re78 ~ treat + age + re75 + educ,
      meth = "nearest",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE,
      tol = 1e-7,
      discard = "both",

      eps = 1e-4, # delta method
      correction_type = "HC0",
      estimand = "ATE",
      treatment = "treat"
    )

    ds_ate <- ds_matching(connections, params_list)
    p <- ds.qq_plot(
      unmatched = "DST",
      matched = "matched_pooled",
      formula = params_list[["matching_form"]],
      datasources = connections
    )

    expect_s3_class(p, "gtable")
  })

test_that(
  "-- eCDF Plotting test --", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <- list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = re78 ~ treat + age + re75 + educ,
      meth = "nearest",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE,
      tol = 1e-7,
      discard = "both",

      eps = 1e-4, # delta method
      correction_type = "HC0",
      estimand = "ATE",
      treatment = "treat"
    )

    ds_ate <- ds_matching(connections, params_list)
    p <- ds.eCDF_plot(
      unmatched = "DST",
      matched = "matched_pooled",
      formula = params_list[["matching_form"]],
      datasources = connections
    )

    expect_s3_class(p, "gtable")
  })


test_that(
  "-- ATE Standard Error test --
  Method: Nearest
  m.order: Smallest
  ratio: 1
  replace: False
  discard: both
  Error: HC0", {

  data_object <- create_test_data(seed = 12345)
  connections <- data_object$connections
  test_data <- data_object$data

  params_list <- list(
    matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
    est_form = re78 ~ treat + age + re75 + educ,
    meth = "nearest",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE,
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC0",
    estimand = "ATE",
    treatment = "treat"
  )

  ds_ate <- ds_matching(connections, params_list)
  avg_non_FL_standard <- avg_local_matching(test_data, params_list)

  expect_equal(avg_non_FL_standard$std.error, ds_ate$se, tolerance=tolerance)
})

test_that(
  "-- ATT Standard Error test --
  Method: Nearest
  m.order: Smallest
  ratio: 1
  replace: False
  discard: both
  Error: HC0", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <- list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = re78 ~ treat + age + re75 + educ,
      meth = "nearest",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE,
      tol = 1e-7,
      discard = "both",

      eps = 1e-4, # delta method
      correction_type = "HC0",
      estimand = "ATT",
      treatment = "treat"
    )

    ds_ate <- ds_matching(connections, params_list)
    avg_non_FL_standard <- avg_local_matching(test_data, params_list)

    expect_equal(avg_non_FL_standard$std.error, ds_ate$se, tolerance=tolerance)
  })
test_that(
  "-- ATC Standard Error test --
  Method: Nearest
  m.order: Smallest
  ratio: 1
  replace: False
  discard: both
  Error: HC0", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <- list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = re78 ~ treat + age + re75 + educ,
      meth = "nearest",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE,
      tol = 1e-7,
      discard = "both",

      eps = 1e-4, # delta method
      correction_type = "HC0",
      estimand = "ATC",
      treatment = "treat"
    )

    ds_ate <- ds_matching(connections, params_list)
    avg_non_FL_standard <- avg_local_matching(test_data, params_list)

    expect_equal(avg_non_FL_standard$std.error, ds_ate$se, tolerance=tolerance)
  })

test_that(
  "-- ATE Standard Error test --
  Method: Nearest
  m.order: Smallest
  ratio: 2
  replace: False
  discard: both
  Error: HC0", {

  data_object <- create_test_data(seed = 12345)
  connections <- data_object$connections
  test_data <- data_object$data

  params_list <- list(
    matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
    est_form = re78 ~ treat + age + re75 + educ,
    meth = "nearest",
    m.order = "smallest",
    ratio = 2,
    replace = FALSE,
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC0",
    estimand = "ATE",
    treatment = "treat"
  )

  ds_ate <- ds_matching(connections, params_list)
  avg_non_FL_standard <- avg_local_matching(test_data, params_list)

  expect_equal(avg_non_FL_standard$std.error, ds_ate$se, tolerance=tolerance)
})

test_that(
  "-- ATE Standard Error test --
  Method: Nearest
  m.order: Smallest
  ratio: 2
  replace: False
  discard: none
  Error: HC0", {

  data_object <- create_test_data(seed = 12345)
  connections <- data_object$connections
  test_data <- data_object$data

  params_list <-   list(
    matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
    est_form = re78 ~ treat + age + re75 + educ,
    meth = "nearest",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "none",

    eps = 1e-4, # delta method
    correction_type = "HC0",
    estimand = "ATE",
    treatment = "treat"
  )

  ds_ate <- ds_matching(connections, params_list)
  avg_non_FL_standard <- avg_local_matching(test_data, params_list)

  expect_equal(avg_non_FL_standard$std.error, ds_ate$se, tolerance=tolerance)
})

test_that(
  "-- ATE Standard Error test --
  Method: Optimal
  m.order: Smallest
  ratio: 1
  replace: False
  discard: both
  Error: HC0", {

  data_object <- create_test_data(seed = 12345)
  connections <- data_object$connections
  test_data <- data_object$data

  params_list <-   list(
    matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
    est_form = re78 ~ treat + age + re75 + educ,
    meth = "optimal",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC0",
    estimand = "ATE",
    treatment = "treat"
  )

  ds_ate <- ds_matching(connections, params_list)
  avg_non_FL_standard <- avg_local_matching(test_data, params_list)

  expect_equal(avg_non_FL_standard$std.error, ds_ate$se, tolerance=tolerance)
})

test_that(
  "-- ATE Standard Error test --
  Method: Optimal
  m.order: Smallest
  ratio: 1
  replace: False
  discard: both
  Error: HC1", {

  data_object <- create_test_data(seed = 12345)
  connections <- data_object$connections
  test_data <- data_object$data

  params_list <-   list(
    matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
    est_form = re78 ~ treat + age + re75 + educ,
    meth = "optimal",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC1",
    estimand = "ATE",
    treatment = "treat"
  )

  ds_ate <- ds_matching(connections, params_list)
  avg_non_FL_standard <- avg_local_matching(test_data, params_list)

  expect_equal(avg_non_FL_standard$std.error, ds_ate$se, tolerance=tolerance)
})

test_that(
  "-- ATE Standard Error test --
  Method: Nearest
  m.order: Smallest
  ratio: 1
  replace: False
  discard: both
  Error: HC2", {

  data_object <- create_test_data(seed = 12345)
  connections <- data_object$connections
  test_data <- data_object$data

  params_list <-   list(
    matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
    est_form = re78 ~ treat + age + re75 + educ,
    meth = "optimal",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC2",
    estimand = "ATE",
    treatment = "treat"
  )

  ds_ate <- ds_matching(connections, params_list)
  avg_non_FL_standard <- avg_local_matching(test_data, params_list)

  expect_equal(avg_non_FL_standard$std.error, ds_ate$se, tolerance=tolerance)
})

test_that(
  "-- ATE Standard Error test --
  Method: Optimal
  m.order: Smallest
  ratio: 2
  replace: False
  discard: both
  Error: HC3", {

  data_object <- create_test_data(seed = 12345)
  connections <- data_object$connections
  test_data <- data_object$data

  params_list <-   list(
    matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
    est_form = re78 ~ treat + age + re75 + educ,
    meth = "optimal",
    m.order = "smallest",
    ratio = 1,
    replace = FALSE, # replace need to not use cluster robustness
    tol = 1e-7,
    discard = "both",

    eps = 1e-4, # delta method
    correction_type = "HC3",
    estimand = "ATE",
    treatment = "treat"
  )

  ds_ate <- ds_matching(connections, params_list)
  avg_non_FL_standard <- avg_local_matching(test_data, params_list)

  expect_equal(avg_non_FL_standard$std.error, ds_ate$se, tolerance=tolerance)
})

test_that(
  "-- Subclass Matching Quantiles Comparison--
  Method: subclass
  n_subclasses: 3
  m.order: Smallest
  ratio: 1
  replace: False
  discard: both
  Estimand: ATE
  Error: HC3", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <-   list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = "re78 ~ treat",
      meth = "subclass",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE, # replace need to not use cluster robustness
      tol = 1e-7,
      discard = "none",

      eps = 1e-4, # delta method
      correction_type = "HC3",
      estimand = "ATE",
      n_subclass = 3,
      treatment = "treat"
    )

    ds_summary <- ds_matching_subclass(
      connections,
      params_list
    )
    non_FL_standard_summary <- local_matching_subclass(test_data, params_list)

    # compare the subclass summaries
    # expect_equal(sum(abs(t(ds_summary$`Sample Sizes`) - non_FL_standard_summary[, 1:(ncol(non_FL_standard_summary)-1)])), 0, tolerance=tolerance)
    # compare qunatiles
    expect_equal(max(ds_summary$quantiles - non_FL_standard_summary$q.cut), 0, tolerance=tolerance)
  })

test_that(
  "-- Subclass Matching Quantiles Comparison--
  Method: subclass
  n_subclasses: 4
  m.order: Smallest
  ratio: 1
  replace: False
  discard: both
  Estimand: ATT
  Error: HC3", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <-   list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = "re78 ~ treat",
      meth = "subclass",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE, # replace need to not use cluster robustness
      tol = 1e-7,
      discard = "none",

      eps = 1e-4, # delta method
      correction_type = "HC3",
      estimand = "ATT",
      n_subclass = 4,
      treatment = "treat"
    )

    ds_summary <- ds_matching_subclass(
      connections,
      params_list
    )
    non_FL_standard_summary <- local_matching_subclass(test_data, params_list)

    # compare the subclass summaries
    # expect_equal(sum(abs(t(ds_summary$`Sample Sizes`) - non_FL_standard_summary[, 1:(ncol(non_FL_standard_summary)-1)])), 0, tolerance=tolerance)
    expect_equal(max(ds_summary$quantiles - non_FL_standard_summary$q.cut), 0, tolerance=tolerance)
  })
test_that(
  "-- Subclass Matching Quantiles Comparison--
  Method: subclass
  n_subclasses: 5
  m.order: Smallest
  ratio: 1
  replace: False
  discard: both
  Estimand: ATT
  Error: HC3", {

    data_object <- create_test_data(seed = 12345)
    connections <- data_object$connections
    test_data <- data_object$data

    params_list <-   list(
      matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
      est_form = "re78 ~ treat",
      meth = "subclass",
      m.order = "smallest",
      ratio = 1,
      replace = FALSE, # replace need to not use cluster robustness
      tol = 1e-7,
      discard = "none",

      eps = 1e-4, # delta method
      correction_type = "HC3",
      estimand = "ATC",
      n_subclass = 5,
      treatment = "treat"
    )

    ds_summary <- ds_matching_subclass(
      connections,
      params_list
    )
    non_FL_standard_summary <- local_matching_subclass(test_data, params_list)

    # compare the subclass summaries
    # expect_equal(sum(abs(t(ds_summary$`Sample Sizes`) - non_FL_standard_summary[, 1:(ncol(non_FL_standard_summary)-1)])), 0, tolerance=tolerance)
    expect_equal(max(ds_summary$quantiles - non_FL_standard_summary$q.cut), 0, tolerance=tolerance)
  })



# data_object <- create_test_data(seed = 12345)
# connections <- data_object$connections
# test_data <- data_object$data
#
# params_list <- list(
#   matching_form = treat ~ age + educ + married + nodegree + re74 + re75,
#   est_form = re78 ~ treat + age + re75 + educ,
#   meth = "nearest",
#   m.order = "smallest",
#   ratio = 1,
#   replace = FALSE,
#   tol = 1e-7,
#   discard = "both",
#
#   eps = 1e-4, # delta method
#   correction_type = "HC0",
#   estimand = "ATE",
#   treatment = "treat"
# )
#
# # normal matching
# norm_match <- MatchIt::matchit(formula = params_list[["matching_form"]],
#                                data = test_data,
#
#                                method = params_list[["meth"]],
#                                m.order = params_list[["m.order"]],
#                                ratio = params_list[["ratio"]],
#                                replace = params_list[["replace"]],
#                                tol = params_list[["tol"]],
#                                discard = params_list[["discard"]])
# # summary of unmatched/matched results
# norm_summary <- summary(norm_match, un = T, improvement = T)
#
# # federated matching
# fed_match <- ds.matchit(form = params_list[["matching_form"]],
#                         data = "DST",
#                         privacy = "norm",
#                         sd = 0,
#                         newobj = "matched_pooled",
#
#                         method = params_list[["meth"]],
#                         m.order = params_list[["m.order"]],
#                         ratio = params_list[["ratio"]],
#                         replace = params_list[["replace"]],
#                         tol = params_list[["tol"]],
#                         discard = params_list[["discard"]],
#                         datasources = connections)
#
# # federated summary
# fed_summary.combined <- ds.match_summary(unmatched_obj = "DST",
#                                          matched_obj = "matched_pooled",
#                                          type = "combined",
#                                          bin_num = 429,
#                                          treatment = params_list[["treatment"]],
#                                          datasources = connections)
#
#
#
#
#
#
# vars_comparison <- all.vars(params_list[["matching_form"]])
# vars_comparison <- vars_comparison[vars_comparison != params_list[["treatment"]]]
# vars_comparison <- c("distance", vars_comparison)
#
# # select rows with same row names
# summed_comparison <- fed_summary.combined$`Summary of Balance for Matched Data:`[vars_comparison, ]
# compared_vec <- summed_comparison$Means.Treated - norm_summary$sum.matched[, "Means Treated"]

gen_data_correct <- function(n=1000, seed=1) {
  set.seed(seed)

  X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
  # induce confounding
  linps <- -1.5 + 0.5*X[,1] + 0.2*X[,2] + 0.2*X[,4] - 0.2*X[,5] - 0.2*X[,7] - 0.2*X[,8]
  A <- rbinom(n, 1, plogis(linps))

  # outcome with constant treatment effect tau
  tau <- 1
  Y <- -2.8 + tau*A + 1.5*X[,1] + 3*X[,2] + 5*X[,3] + 4*X[,4] + 4*X[,5] + 4*X[,6] + rnorm(n,0,5)

  data.frame(A, X, Y)
}

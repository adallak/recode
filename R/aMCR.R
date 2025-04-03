#' Implements sparse multivariate regression
#'
#' This function implements joint estimation of sparse multivariate
#' regression and conditional graphical models as in Wang (2015)
#' (https://www3.stat.sinica.edu.tw/statistica/oldpdf/A25n31.pdf)
#'
#' @param Y  multivariate response
#' @param X  data matrix
#' @param lam.vec vector of lambda values. None by default.
#' @param nfolds number of folds; 10 by default
#' @param penalty penalty type; lasso, MCP and SCAF. lasso is the default
#' @param symm.rule the rule to make Gamma_hat matrix symmetric, if True `or` rule is used, and `and` rule otherwise
#' @param seed   random seed
#'
#' @return
#' *`B_hat` - Estimated coefficient matrix
#' *`Gamma_hat` - Estimated gamma matrix with the same sparsity as the precision matrix
#' *`B_0_hat` - Estimated initial coefficient matrix
#' *`lambda.min` - the vector of lambda values selected from cross-validation
#' @export
aMCR <- function(Y, X, lam.vec = NULL, nfolds = 10,
                penalty = c("lasso", "MCP", "SCAD"),
                symm.rule = TRUE , seed = 123){

  penalty = match.arg(penalty)
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  if (n != dim(Y)[1]){
    stop("Y and X should have the same number of observations")
  }
  B_0_hat = matrix(0, nrow = p, ncol = q)
  ## Compute initial matrix B_0
  for (i in 1:q){
      cv.lasso = glmnet::cv.glmnet(x = X, y = Y[,i], lambda = lam.vec)
      B_0_hat[,i] = as.matrix(coef(cv.lasso))[2:(p+1)]
  }
  if (is.null(lam.vec)){
      lam.vec = cv.lasso$lambda
  }
  ## run aMCR
  B_hat = matrix(0, nrow = p, ncol = q)
  Gamma_hat = matrix(0, nrow = q, ncol = q)
  lambda.min = c()
  results = c()
  for (k in 1:q){
      mat = Y[,-k] - X %*% B_0_hat[,-k]
      X_aug = cbind(X, mat)
      if (penalty == "lasso"){
        result = glmnet::cv.glmnet(x = X_aug, y = Y[,k], lambda = lam.vec)
      } else {
        result =  ncvreg::cv.ncvreg(X_train, y_train,
                                    penalty = penalty, seed = seed)
      }
      coef = as.matrix(coef(result))[2,(p + q + 1)]
      B_hat[,k] = coef[1:p]
      Gamma_hat[,k] = coef[(p+1):q]
      lambda.min[i] = result$lambda.min
      results[i] = result
  }
  if (isTRUE(symm.rule)){
      mask =check_zero_entries(Gamma_hat)
  }
  else {
      mask =cons_check_zero_entries(Gamma_hat)
  }
  Gamma_hat = mask * Gamma_hat
  return(list("B_hat" = B_hat, "Gamma_hat" = Gamma_hat,
              "B_0_hat" = B_0_hat,
              "lambda.min" = lambda.min, "lambda" = lam.vec,
              "model" = "aMCR"))
}

check_zero_entries <- function(mat) {
    return((mat == 0) | (t(mat) == 0))
}

cons_check_zero_entries <- function(mat) {
  return((mat == 0) & (t(mat) == 0))
}




compute_var_ij <- function(x, i, j) {
  # X: n x p matrix (n observations, p variables)
  n <- nrow(x)
  p <- ncol(x)

#  term1 = sum((x[,i] - mean(x[,i]))^2 * (x[,j] - mean(x[,j]))^2)
#  term2 = (sum((x[,i] - mean(x[,i])) *(x[,j] - mean(x[,j]))))^2 / n

  return(n / (n-1)^3 * sum(((x[,i] - mean(x[,i])) * (x[,j] - mean(x[,j])) -
                             1/n * sum((x[,i] - mean(x[,i])) * (x[,j] - mean(x[,j]))))^2))
}

compute_cov_ij_lm <- function(x, i, j, l, m){
  n <- nrow(x)
  p <- ncol(x)

  term__ij = ((x[,i] - mean(x[,i])) * (x[,j] - mean(x[,j])) -
      1/n * sum((x[,i] - mean(x[,i])) * (x[,j] - mean(x[,j]))))
  term__lm = ((x[,l] - mean(x[,l])) * (x[,m] - mean(x[,m])) -
      1/n * sum((x[,l] - mean(x[,l])) * (x[,m] - mean(x[,m]))))

    return(n / (n-1)^3 * sum(term__ij * term__lm))
}

s_ij <- function(x, i, j){
  n <- nrow(x)
  p <- ncol(x)
  return(n/(n-1) *  1/n * sum((x[,i] - mean(x[,i])) * (x[,j] - mean(x[,j]))))
}

compute_f_ij <- function(x, i, j) {
  # X: n x p matrix (n observations, p variables)
  return(1/2 * (sqrt(s_ij(x,j,j) / s_ij(x, i, i)) *
                  compute_cov_ij_lm(x, i, i, i, j) +
                sqrt(s_ij(x,i,i) / s_ij(x, j, j)) *
                  compute_cov_ij_lm(x, j, j, i, j))
         )
}


numerator <- function(x){
  # X: n x p matrix (n observations, p variables)
  n <- nrow(x)
  p <- ncol(x)
  sum = 0
  r_bar = sum((cor(x) - diag(p)))/ (p * (p-1))
  for (i in 1:p){
    for (j in 1:p){
      if (i != j){
        sum = sum + compute_var_ij(x, i, j) - r_bar * compute_f_ij(x, i, j)
      }
    }
  }
  return(sum)
}

denominator <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  sum = 0
  r_bar = sum((cor(x) - diag(p)))/ (p * (p-1))
  for (i in 1:p){
    for (j in 1:p){
      if(i != j){
        sum = sum + (s_ij(x,i, j) - r_bar *  sqrt(s_ij(x,i,i) * s_ij(x,j,j)))^2
      }
    }
  }
  return(sum)
}

estimate_lambda_F <- function(x){
  return(max(0, min(1, numerator(x) / denominator(x))))
}

#' Implements shrinkage estimation of the covariance matrix using case F as in
#' Schaffer and Strimmer (2005)
#'
#' This function implements shrinkave covariance matrix estimation
#' using case F as in Schaffer and Strimmer (2005) (https://strimmerlab.github.io/publications/journals/shrinkcov2005.pdf)
#'
#' @param x  data matrix
#'
#' @return
#' *`Sigma` - shrinkage covariance matrix
#' @export
cov_shrink_F <- function(x){
  n <- nrow(x)
  p <- ncol(x)

  res = corpcor::var.shrink(x)
  D = diag(sqrt(res))

  lambda_F   = estimate_lambda_F(x)
  cat("Estimating optimal shrinkage intensity lambda for case F :", lambda_F, "\n")
  #cov_X = cov(S)
  r_bar = sum((cor(x) - diag(p)))/ (p * (p-1))
  ##rho = (sum(cov_X) - sum(diag(cov_X))) / (p * (p - 1))
  mat_T = r_bar * diag(p) + (1 - r_bar) * matrix(1, p, p)
  R = lambda_F * mat_T + (1 - lambda_F) * cor(x)

  sigma = D %*% R %*% D
  attr(sigma, "lambda.var") = attr(res, "lambda.var")
  attr(sigma, "lambda") = lambda_F

  return (sigma)
}



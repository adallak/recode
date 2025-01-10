#' Generate data from latent variable graphical model
#'
#'
#' @param p_o  number of observed variables; default is 10
#' @param p_l  number of latent variables; default is 5
#' @param N    number of observations; default is 150
#' @param prob_o probability two observed variables are connected; default is 0.4
#' @param prob_l probability two latent variables are connected; default is 0.2
#' @param prob_ol probability the latent and observed variables are connected; default is 0.2
#' @param lwb  lower value for the elements of the inverse covariance matrix
#' @param seed random seed
#' @returns
#' * `X` N x p_o data matrix
#' * `theta_o` true observed inverse covariance matrix
#' * `theta_obs` observed inverse covariance matrix
#' @export
genlatentdata <- function(p_o = 10, p_l = 5, N = 150,
                          prob_o = 0.4, prob_l = 0.2, prob_ol = 0.7,
                          lwb = 0.4, seed = 1234) {
   set.seed(seed)
   sparse = rbinom((p_o) * (p_o), 1, prob_o) #*  runif(p_o * p_o, lwb, 1)
   theta_o = matrix(sparse, nrow = p_o, ncol = p_o)
   theta_o = theta_o + t(theta_o)
   theta_o[theta_o>1] = 1
   diag(theta_o) = 1
   ee = min(Re(eigen(theta_o, only.values = T)$values))
   diag(theta_o) = 1 + ifelse(ee < 0, -ee + 0.3, 0.3)
   if (p_o <= 0) {
      stop("number of observed variables should be positive")
   }
   if (p_l < 0) {
      stop("number of latent variables should be positive")
   }
   if (N <= 0) {
      stop("number of observations should be positive")
   }

   if (p_l > 0) {
      sparse = rbinom((p_l) * (p_l), 1, prob_l) #*  runif(p_l *p_l, lwb, 1)
      theta_l = matrix(sparse, nrow = p_l, ncol = p_l)
      theta_l = theta_l + t(theta_l)
      diag(theta_l) = 1
      theta_l[theta_l > 1] = 1
      ee = min(Re(eigen(theta_l, only.values = T)$values))
      diag(theta_l) = 1 + ifelse(ee < 0, -ee + 0.3, 0.3)
      sparse = rbinom((p_o) * (p_l), 1, prob_ol) * 0.2
      theta_ol = matrix(sparse, nrow = p_o, ncol = p_l)
      theta_lo = t(theta_ol)
      theta_cond_o =  theta_o - theta_ol %*% solve(theta_l) %*% theta_lo
      #ee = min(Re(eigen(theta_cond_o, only.values = T)$values))
      #diag(theta_cond_o) = diag(theta_cond_o) + ifelse(ee < 0, -ee + 0.3, 0.3)
      sigma_cond_o = solve(theta_cond_o)
      mu = matrix(0, nrow = p_o)
      X =  matrix(MASS::mvrnorm(N, mu, sigma_cond_o), nrow =N,
                  ncol = p_o)
   }
   else {
      theta_cond_o =  NULL
      sigma_o = solve(theta_o)
      mu = matrix(0, nrow = p_o)
      X =  matrix(MASS::mvrnorm(N, mu, sigma_o), nrow =N,
                  ncol = p_o)
   }
   return(list(X = X, theta_o = theta_o, theta_obs = theta_cond_o ))
}


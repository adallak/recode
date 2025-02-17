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


#' Generate data from AD or Star shape undirected graph
#'
#'
#' @param n    number of observations
#' @param p    number of variables
#' @param type type of the undirected graph
#' @param variance variance of the generated data; default is 1
#' @param rho  defines edge weights
#' @param AD_order order of the AD process; ignored for `type(star)`
#' @param seed random seed
#' @returns
#' *`omega` p x p precision matrix
#' * `data` n x p data

#' @export
generatedata <- function(n, p, type = c("star", "AD1"),
                         variance = 1, rho = 0.5, AD_order = 1,
                         seed = 123) {
   # n: Number of observations
   # p: Number of variables (dimensions)
   # central_variance: Variance of the central node
   # edge_weight: Precision weight connecting the central node to others
   set.seed(seed)
   if (AD_order >= p) {
      stop("Error: 'AD_order' must be less than 'p' (number of variables).")
   }
   if (n <= 0) {
      stop("Error: 'n' must be positive number.")
   }
   if (p <= 0) {
      stop("Error: 'p' must be positive number.")
   }

   type = match.arg(type)
   # Initialize precision matrix
   Omega <- matrix(0, p, p)

   if(type == "star") {
      # Set central node index (1st variable)
      Omega[1, 1] <- central_variance + (p - 1) * edge_weight
      for (i in 2:p) {
         Omega[i, i] <- 1  # Variance of peripheral nodes
         Omega[i, 1] <- -edge_weight  # Connection to the central node
         Omega[1, i] <- -edge_weight  # Symmetric connection
      }
   }
   else {
      for (i in 1:p) {
         for (j in max(1, i - AD_order): (i - 1)) {  # Only connect up to t_order previous vars
            Omega[i, j] <- -rho^(i - j)  # Decaying correlation
            Omega[j, i] <- -rho^(i - j)  # Symmetric assignment
         }
         Omega[i, i] <- 1 + sum(rho^(1:min(AD_order, i - 1)))  # Adjust diagonal elements
      }

   }
   # Compute covariance matrix
   Sigma <- solve(Omega)

   # Generate data from multivariate normal distribution
   data <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

   return(list("Omega" = Omega, "data" = data))
}


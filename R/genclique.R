#' Generate a block diagonal covariance matrix
#' This function is borrowed from `spcov` package.
#' For details, see https://cran.r-project.org/web/packages/spcov/index.html.
#'
#'
#' @param ncliques  number of blocks; default is 5
#' @param cliquesize  size of each block; default is 5
#' @param theta    magnitude of non-zeros elements; default is 1
#' @param seed     random seed
#' @returns
#' * `Sigma` p x p covariance matrix
#' * `A` p x p square root of  `Sigma`
#' @export
genclique <- function (ncliques = 5, cliquesize = 5, theta = 1, seed = 1234)
{
   set.seed(seed)
   if ((ncliques <= 0) | (ncliques %% 1 != 0)) {
      stop("number of blocks should be positive integer")
   }
   if ((cliquesize <= 0) | (cliquesize %% 1 != 0)) {
      stop("size of each block should be positive integer")
   }
   if (theta <= 0) {
      stop("magnitude of non-zeros elements should be positive")
   }
   p <- ncliques * cliquesize
   sizes <- rep(cliquesize, ncliques)
   Sigma <- matrix(0, p, p)
   lims <- c(0, cumsum(sizes))
   for (i in seq(ncliques)) {
      ii <- (lims[i] + 1):lims[i + 1]
      signs <- 2 * (stats::rnorm(sizes[i] * (sizes[i] - 1)/2) <
                       0) - 1
      Sigma[ii, ii][upper.tri(Sigma[ii, ii])] <- signs * theta
   }
   Sigma <- Sigma + t(Sigma)
   eig <- eigen(Sigma, symmetric = T)
   shift <- (max(eig$val) - p * min(eig$val))/(p - 1)
   cat("Shifting eigenvalues by ", shift, fill = T)
   diag(Sigma) <- diag(Sigma) + shift
   A <- eig$vect %*% diag(sqrt(eig$val + shift)) %*% t(eig$vect)
   return(list(Sigma = Sigma, A = A, shift = shift))
}


#' Generate data from a block diagonal covariance matrix
#'
#' @param p  number of features
#' @param n  number of observations
#' @param ncliques  number of blocks; default is 5
#' @param theta    magnitude of non-zeros elements; default is 1
#' @param seed     random seed
#' @returns
#' * `X` n x p data matrix
#' * `S` p x p sample covariance matrix
#' * `Sigma` p x p true covariance matrix
#' * `A` p x p square root of  `Sigma`
#' @export
gendataclique <- function (p = 20, n = 100, ncliques = 5,
                       theta = 1, seed = 1234)
{
   set.seed(seed)
   cliquesize = p / ncliques
   if ((ncliques <= 0) | (ncliques %% 1 != 0)) {
      stop("number of blocks should be positive integer")
   }
   if ((p <= 0) | (p %% 1 != 0)) {
      stop("number of features be positive integer")
   }
   if ((n <= 0) | (n %% 1 != 0)) {
      stop("number of observations be positive integer")
   }
   if (theta <= 0) {
      stop("magnitude of non-zeros elements should be positive")
   }
   genc = genclique(ncliques, cliquesize, theta, seed)
   Sigma = genc$Sigma
   A     = genc$A
   shift = genc$shift
   X <- matrix(stats::rnorm(n * p), ncol=p) %*% genc$A
   S <- stats::var(X)
   return(list(X = X, S = S, Sigma = Sigma, A = A, shift = shift))
}

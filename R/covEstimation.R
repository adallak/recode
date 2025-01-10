#' Implement covariance matrix estimation
#'
#' This function implements covariance matrix estimation
#' using:
#' - Thresholding method proposed in Cai and Liu (2011) and Bickel and Levina (2008b)
#' - Banding method proposed in Bickel and Levina (2008a)
#' - Tapering method proposed in Cai et al (2010)
#'
#' @param X  data matrix
#' @param k  tuning parameter
#' @param method if `method = "thresholding` and the `thr.method` is `adaptive` implements
#'          method by Cai and Liu (2011), otherwise method by Bickel and Levina (2008b).
#' @param thr.method method for thresholding
#' @param demean demeanize the data. True by default
#' @param scale standardize the data. False by defaults.
#'
#' @return
#' *`Sigma` - Estimated covariance matrix
#' *`S`     - Sample covariance matrix
#' @export
covEstimation <- function(X, k, method = c("thresholding", "banding","tapering"),
                          thr.method = c("adaptive", "soft", "hard"),
                          demean = FALSE, scale = FALSE){
      method = match.arg(method)
      thr.method = match.arg(thr.method)
      S = cov(scale(X,center = demean, scale = scale))
      p = ncol(X)
      n = nrow(X)
      if (method == "banding"){
            Sigma = threshmat(S, k)
      }
      if (method == "tapering"){
            Sigma = S * tapermat(S, k)
      }
      if (method == "thresholding"){
            if (thr.method == "adaptive")
            {
                  diag_S = diag(S)
                  theta = variances_cov_matrix(X, S)
                  Lambda = k * sqrt(theta * log(p)/n)
                  Sigma = pmax(abs(S)-Lambda, 0)*sign(S)
                  diag(Sigma) = diag_S
            }
            else if (thr.method == "soft")
            {
                  diag_S = diag(S)
                  Sigma = pmax(abs(S)-k, 0)*sign(S)
                  diag(Sigma) = diag_S
            }
            else{
                  bool = abs(S) > k
                  Sigma = bool * S
                  diag(Sigma) = diag(S)
            }
      }

      return(list("Sigma" = Sigma, "S" = S))
}

threshmat <- function(M, k){
      p = nrow(M)
      for (i in 1:p){
            for (j in 1:p){
                  if (abs(i - j) > k){
                        M[i,j] = 0
                  }
            }
      }
      return(M)
}

tapermat <- function(M, k){
      p = nrow(M)
      W = matrix(0, p,p)
      for (i in 1:p){
            for (j in 1:p){
                  if (abs(i - j) <= k){
                        W[i,j] = 2 - 2 * abs(i - j) / k
                  }
                  if (abs(i - j) <= k/2) {
                        W[i,j] = 1
                  }
            }
      }
      return(W)
}

variances_cov_matrix <- function(X, S){
      variances_cov_matrix <- matrix(0, nrow = ncol(S),
                                     ncol = ncol(S))
      n = nrow(X)
      for (i in 1:ncol(S)) {
            for (j in 1:ncol(S)) {
                  cov_entry <- S[i, j]
                  col_mean = colMeans(X)
                  mean_entry <- (X[, i] - col_mean[i]) * (X[, j] - col_mean[j])
                  variances_cov_matrix[i, j] <- sum((cov_entry - mean_entry)^2
                                                    ) / (n)
            }
      }
      return(variances_cov_matrix)
}

#' Implement covariance matrix estimation
#'
#' This function implements covariance matrix estimation
#' using:
#' - Thresholding method proposed in Cai and Liu (2011) and Bickel and Levina (2008b)
#' - Banding method proposed in Bickel and Levina (2008a)
#' - Tapering method proposed in Cai et al (2010)
#'
#' @param X  data matrix
#' @param k.list  list of tuning parameter; default is NULL
#' @param n.k     number of tuning parameters; default is 30
#' @param k.min    minimum value for tuning parameter; default is 0.01
#' @param k.max    maximum value for tuning parameter; no maximum is specified by default
#' @param nfolds   number of folds for CV; default is 5
#' @param method if `method = "thresholding` and the `thr.method` is `adaptive` implements
#'          method by Cai and Liu (2011), otherwise method by Bickel and Levina (2008b).
#' @param thr.method method for thresholding
#' @param demean demeanize the data. True by default
#' @param scale standardize the data. False by defaults.
#'
#' @return
#' *`err_fit` - fitted errors
#' *`cv,k` - selected tuning parameter
#' *`Sigma_fit` - Estimated covariance matrix
#' *`S`     - Sample covariance matrix
#' @export
cv.covEstimation <- function(X, k.list = NULL, n.k = 30,
                             k.min = 0.01, nfolds = 5,
                 method = c("thresholding", "banding", "tapering"),
                 thr.method = c("adaptive", "soft", "hard"),
                 demean = FALSE, scale = FALSE, k.max = NULL) {
      n <- nrow(X)
      p <- ncol(X)
      S = cov(scale(X,center = demean, scale = scale))
      folds = makefolds(n, nfolds)
      method = match.arg(method)
      thr.method = match.arg(thr.method)
      if (is.null(k.list)){
            if (method == "thresholding")
            {
                  if (is.null(k.max))
                  {
                        k.max = k_max(S)
                  }
                  else{
                        k.max = k.max
                  }
                  k.list = gen.list(n.k, k.max, k.min, S)
            }
      }
      else{
                  k.list = seq(1:p)
                  n.k = length(k.list)
            }
      errs_fit <- matrix(NA, n.k, nfolds)
      est_cov = array(NA, dim = c(p,p,n.k))
      for (i in 1:nfolds){
            x_tr <- X[-folds[[i]],]
            meanx <- colMeans(x_tr)
            x_tr <- scale(x_tr, center = meanx, scale = scale)
            S_tr <- cov(x_tr)
            iter = 1
            for (k in k.list){
                 est_cov[,,iter] = covEstimation(x_tr, k, method = method,
                                                 thr.method = thr.method,
                                            demean = demean, scale = scale)$Sigma
                 iter = iter + 1
            }
            x_te <- X[folds[[i]], ]
            x_te <- scale(x_te, center = meanx, scale = scale)
            S_te <- cov(x_te)

            for (j in 1:n.k) {
                  errs_fit[j, i] <- est.error(est_cov[,,j], S_te)
            }
      }
      ibest_fit <- which.min( rowMeans(errs_fit))

      sigma_fit = covEstimation(X, k.list[ibest_fit], method = method, thr.method = thr.method,
                              demean = demean, scale = scale)$Sigma
      return(list(errs_fit = errs_fit, folds = folds, k.list = k.list,
                  cv.k = k.list[ibest_fit], Sigma_fit = sigma_fit, S = S))
}

gen.list <- function(nlam, lam_max, flmin, S){
      ## This function is borrowed from varband package
      lamlist_lin <- lam_max * exp(seq(0, log(flmin), length = nlam/2))
      lamlist_exp <- seq(lam_max - 1e-8, lam_max*flmin - 1e-8, length.out = nlam/2)
      return(sort(unique(c(lamlist_lin, lamlist_exp)), decreasing = T))
}

k_max <- function(S){
      #### This function is borrowed from varband package
      p <- ncol(S)
      sighat <- rep(NA, p-1)
      for (r in seq(2, p)){
            sighat[r-1] <- max(abs(S[(1:(r-1)), r]))/sqrt(S[r, r])
      }
      return(2 * max(sighat))
}

makefolds <- function(n, nfolds) {
      nn <- round(n / nfolds)
      sizes <- rep(nn, nfolds)
      sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
      b <- c(0, cumsum(sizes))
      ii <- sample(n)
      folds <- list()
      for (i in seq(nfolds)) {
            folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
      }
      return(folds)
}

est.error <- function(S.hat, S){
      return(sum((S.hat - S)^2))
}

#' Plot the covariance matrix. Borrowed from `varband` package.
#'
#'
#' @param Mat  covariance matrix
#' @param main  title of the plot
#'
#' @export
plot_mat <- function(Mat, main = NULL){
      #### This function is borrowed from varband package
      tmppar <- par(pty = "s")
      image(sign(t(apply(Mat, 2, rev))), axes = FALSE, col = c("gray50","white","black"), main = main)
      par(tmppar)
}

#' Finds the threshold parameter `k.min` that makes the thresholded covariance
#' matrix positive definite
#'
#'
#' @param X  data matrix
#' @param k.grid  grid of tuning parameter
#' @param thr.method method for thresholding
#' @param demean demeanize the data. True by default
#' @param scale standardize the data. False by defaults.
#'
#' @return
#' *`k.min` - selected threshold
#' *`eigen.store`     - stored minimum eigenvalues for each k in the `k.grid`
#' *`k.grid`          - grid of tuning parameters
#' @export
find.thr.min <- function(X, k.grid, thr.method = c("adaptive", "soft", "hard"),
                          demean = FALSE, scale = FALSE){
   thr.method = match.arg(thr.method)
   i = 1
   eigen.store = c()
   for (k in k.grid){
      est.cov = covEstimation(df, k = k, method = "thresholding",
                              thr.method = thr.method)$Sigma
      eigen.store[i] = min(eigen(est.cov)$values)
      i = i + 1
   }
   min.index = which(eigen.store>0)[1]
   if (is.na(min.index)){
      stop("no positive eigenvalue is found. Try different grid.")
   }

   return(list("k.min" = k.grid[min.index], "eigen.store" = eigen.store, "k.grid" = k.grid))
}

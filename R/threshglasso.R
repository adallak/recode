#' Implement and select threshold for thresholded GLasso
#'
#' This function implements tuning of thresholded GLasso
#' using eBIC, AIC, and BIC
#'
#' @param X  data matrix
#' @param lambda  penalization parameter; default is 0.01
#' @param ngrid number of grids for tuning threshold; default is 30
#' @param min.thresh minimum threshold value; default is 0.01
#' @param max.thresh maximum threshold value; default is maximum non-diagonal value of sample covariance matrix
#' @param grid user specified grid
#' @param tol tolerance for GLasso
#' @param gamma parameter for eBIC; default is 0.5
#' @param diagonal pernalize diagonal
#' @param crit criteria for threshold selection
#' @param maxit number of iterations
#' @param start `warm` or `cold` start to improve convergence
#' @param ...   other options
#'
#' @return estimated (inverse) covariance matrix
#' @export
threshglasso <- function(X = NULL, lambda = 0.01, ngrid = 30, min.thresh = 0.01,
                         max.thresh = NULL, grid = NULL, tol = 1e-5, gamma = 0.5,
                         diagonal = FALSE, crit = c("BIC", "AIC", "eBIC"), maxit = 1000,
                         start= c("warm", "cold"), ...) {
   if (!all(c(ngrid,lambda) > 0)) {
      stop("one of ngrid or lambda is negative!")
   }
   if (!all(c(min.thresh,max.thresh) > 0)) {
      stop("one of min.thresh or max.thresh is negative!")
   }
   p = ncol(X)
   crit = match.arg(crit)
   start = match.arg(start)

   S = (nrow(X) - 1)/nrow(X) * stats::cov(X)
   S_nd = S
   diag(S_nd) = 0
   if (is.null(max.thresh)) {
      max.thresh = max(abs(S_nd))
   }
   if (is.null(grid)) {
      grid = 10^seq(log10(min.thresh), log10(max.thresh), length = ngrid)
   } else {
      if ((any(grid) < 0) | length(grid) <= 1) {
         stop("grid should have length < 1 and contain positive values!")
      }
      grid = sort(grid)
   }
   init = S
   initOmega = diag(p)
   errors = array(0, length(ngrid))
   # update diagonal elements of init, if necessary
   if (diagonal) {
      diag(init) = diag(S) + lambda
   }

   # compute the penalized likelihood precision matrix
   # estimator
   GLASSO = glasso::glasso(s = S, rho = lambda, thr = tol,
                           maxit = maxit, penalize.diagonal = diagonal,
                           start = "warm", w.init = init, wi.init = initOmega,
                           trace = FALSE, ...)


   for (i in 1:length(grid)) {

      # set temporary threshold
      thr_ = grid[i]

      thr_omega = GLASSO$wi
      thr_omega[thr_omega < thr_] = 0
      diag(thr_omega) = diag(GLASSO$wi)

      # compute the observed negative validation loglikelihood
      # (close enoug)
      errors[i] = (nrow(X)) * (sum(thr_omega *
                                      S) - determinant(thr_omega, logarithm = TRUE)$modulus[1])

      # update for crit.cv, if necessary
      if (crit == "AIC") {
         errors[i] = errors[i] + sum(thr_omega != 0)
      }
      if (crit == "BIC") {
         errors[i] = errors[i] + sum(thr_omega != 0) * log(nrow(X))

      }
      if (crit == "eBIC") {
         errors[i] = errors[i] + (sum(thr_omega != 0) * log(nrow(X)) +
                                     4 * sum(thr_omega != 0) * gamma * log(ncol(X)))


      }
   }

   # determine optimal tuning parameters
   best_thr = grid[which.min(errors)]
   error = errors[which.min(errors)]
   thr_glasso = GLASSO$wi
   thr_glasso[thr_glasso < best_thr] = 0

   return(list(thr_glasso = thr_glasso, thr = best_thr, min.error = error,
               error = errors, grid = grid))
}


#' Implement and select threshold for thresholded Lasso
#'
#' This function implements tuning of thresholded Lasso
#' using cross-validation
#'
#' @param x  data matrix
#' @param y dependent variable
#' @param lambda  penalization parameter; default is 0.01
#' @param ngrid number of grids for tuning threshold; default is 30
#' @param min.thresh minimum threshold value; default is 0.01
#' @param max.thresh maximum threshold value; default is maximum non-diagonal value of sample covariance matrix
#' @param grid user specified grid
#' @param nfold number of folds for CV
#' @param seed random seed
#' @param fold_ids  (optional) vector of length n specifiying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored
#' @param ...   other glmnet options
#'
#' @return coefficient matrix
#' @export
cv.threshlasso <- function(x, y, lambda = NULL, ngrid = 30, min.thresh = NULL,
                           max.thresh = NULL, grid = NULL, nfold = 5, seed = 1234,
                           fold_ids = NULL, ...) {
   set.seed(seed)
   if (!all(c(ngrid,lambda) > 0)) {
      stop("one of ngrid or lambda is negative!")
   }
   if (!all(c(min.thresh,max.thresh) > 0)) {
      stop("one of min.thresh or max.thresh is negative!")
   }
   if (!is.null(colnames(x))){
      col.names = colnames(x)
   }
   x = as.matrix(x)
   colnames(x) = col.names
   p = ncol(x)
   n = nrow(x)

   if (is.null(max.thresh)) {
      Ymean <- mean(y)
      Ytilde <- y - Ymean
      Xmeans <- colMeans(x)
      Xprime <- x - matrix(colMeans(x), n, ncol(x), byrow = T)
      # Scale X
      weights <- sqrt(colSums(Xprime^2) / n)
      Xtilde <- Xprime %*% diag(1 / weights)
      max.thresh <- max(abs(crossprod(Xtilde, Ytilde)))/n
   }
   if(is.null(lambda)) {
      lambda = glmnet::cv.glmnet(x, y)$lambda.min
   } else {
      lambda = lambda
   }

   if(is.null(min.thresh)) {
      min.thresh = lambda
   }

   if (is.null(grid)) {
      grid = 10^seq(log10(min.thresh), log10(max.thresh), length = ngrid)
   } else {
      if ((any(grid) < 0) | length(grid) <= 1) {
         stop("grid should have length < 1 and contain positive values!")
      }
      grid = sort(grid)
   }
   n_grid = length(grid)
   if (is.null(fold_ids)){
      # Sample splits
      fold_id_n <- sample(rep(seq_len(nfold), length.out = n))
   }else{
      fold_id_n <- fold_ids
      nfold <- max(fold_ids)
   }
   # Matrix to store the residuals
   rCV <- matrix(0, n, n_grid)
   # Matrix to store error means per fold
   CVfolds <- matrix(0, nfold, n_grid)
   # Calculate Lasso for each fold removed
   for (i in 1:nfold){
      # Apply lasso with ith fold removed
      out_fold <- fitLASSO_seq(X = x[fold_id_n != i,], y = y[fold_id_n != i],
                               lambda= lambda, grid = grid, ...)
      beta_mat <- out_fold$beta_mat
      rownames(beta_mat) = col.names
      beta0_vec <- out_fold$beta_0
      # Calculate rCV - residuals for ith fold
      rCV[fold_id_n == i, ] = matrix(y[fold_id_n == i], nrow = sum(fold_id_n == i),
                                     ncol = n_grid) -
         matrix(beta0_vec, nrow = sum(fold_id_n == i),
                ncol = n_grid, byrow = T) -  x[fold_id_n == i, ]%*%beta_mat
      # Calculate CVfolds - means of nfold fold
      CVfolds[i, ] = colMeans(rCV[fold_id_n == i, ]^2)
   }
   # Calculate CV mean and SE
   # Square all the calculated residuals
   rCV <- rCV^2
   # Calculate the CV(lambda)
   cvm <- colMeans(rCV)
   # Calculate the SE(CV(lambda))
   cvse <- apply(CVfolds, 2, sd)/sqrt(nfold)
   # Find threshold_min
   id <- which.min(cvm)
   thresh_min <- grid[id]
   # Find thresh_1SE
   thresh_1se <- max(grid[cvm <= cvm[id] + cvse[id]])
   ## Estimate final
   out <- glmnet::glmnet(x= x, y = y, lambda = lambda, ...)
   beta.coef = coef(out, lambda = lambda)
   beta.lasso.vec <- beta.lasso[-1,]  # Remove intercept
   active_set <- which(abs(beta.lasso.vec) > thresh_min)
   if(length(active_set) > 0) {
     # Refit with OLS on selected variables
     X_active <- x[, active_set, drop = FALSE]
     lm_fit <- lm(y ~ X_active)

     # Initialize beta vector with zeros
     beta.thr <- rep(0, p)
     names(beta.thr) <- colnames(x)

     # Fill in the refitted coefficients
     beta.thr[active_set] <- coef(lm_fit)[-1]  # Exclude intercept
     beta0.thr <- coef(lm_fit)[1]  # Intercept
   } else {
     # No variables selected, return all zeros
     beta.thr <- rep(0, p)
     names(beta.thr) <- colnames(x)
     beta0.thr <- mean(y)  # Just the mean of y
   }

   return(list(grid = grid,
               beta0.thr = beta0.thr, beta.thr = beta.thr,
               beta_mat = beta_mat,
               beta0_vec = beta0_vec, lambda = lambda,
               fold_ids = fold_id_n, thresh_min = thresh_min,
               thresh_1se = thresh_1se, cvm = cvm, cvse = cvse))
}


#' Predicts thresholded Lasso
#'
#'
#' @param thresh  object from the `cv.threshlasso`
#' @param test_x  new dataset for prediction
#'
#' @return vector of predictions
#' @export
predict_threshlasso <- function(thresh, test_x) {
  # Convert to matrix if necessary
  if (!is.matrix(test_x)) {
    test_x <- as.matrix(test_x)
  }

  # Validate dimensions
  if (ncol(test_x) != length(thresh$beta.thr)) {
    stop(paste("test_x has", ncol(test_x), "columns but model expects",
               length(thresh$beta.thr)))
  }

  # Get coefficients
  beta_0 <- thresh$beta0.thr
  beta_vec <- thresh$beta.thr  # More accurate name

  # Make predictions
  predictions <- beta_0 + test_x %*% beta_vec

  # Return as vector (drop matrix dimension)
  return(as.vector(predictions))
}


fitLASSO_seq <- function(X, y, lambda, grid, ...){
   # Calculate solution for each lambda
   n.thresh = length(grid)
   p = dim(X)[2]
   beta_mat <- matrix(0, p, n.thresh)
   beta_0 <- matrix(0, 1, n.thresh)
   out <- glmnet::glmnet(x= X, y = y, lambda = lambda, ...)
   for (i in 1:n.thresh){
      beta.coef = coef(out, s = lambda)
      beta_0[,i] = beta.coef[1,]
      beta.coef = beta.coef[-1,]
      S = which(abs(beta.coef) > grid[i])
      zero = which(abs(beta.coef) <= grid[i])
      beta.coef[zero] = 0
      if(length(S) > 0){
         beta = solve(t(X[,S]) %*% X[, S], t(X[,S]) %*% y)
         beta.coef[S] = beta
      }
         beta_mat[ , i] <- beta.coef
   }
   return(list("beta_0" = beta_0, "beta_mat" =  beta_mat))
}

threshlasso <- function(X, y, twostep = FALSE, seed = 123, ...) {
   if(nrow(X) != length(y)){
      stop("X and y should have the same number of dimensions")
   }
   set.seed(seed)
   lasso_init = glmnet::cv.glmnet(X, y)
   lambda.init = lasso_init$lambda.min
   beta.coef = coef(lasso_init, s = lambda.init)
   beta = beta.coef
   S_0 = which(abs(beta.coef) > 4 * lambda.init)
   if (isFALSE(twostep)){
      zero = which(beta.coef < 4 * lambda.init)
      beta[zero] = 0
      if(length(S_0) > 0){
         beta_1_I = solve(t(X[,S_0]) %*% X[, S_0], t(X[,S_0]) %*% y)
         beta[S_0] = beta_1_I
      }
   } else{
      for (i in 0:1){
         t = 4 * lambda.init * sqrt(length(S_0))
         S = which(beta.coef >= t)
         zero = which(beta.coef < t)
         beta[zero] = 0
         if(length(S) > 0){
            beta_1_I = solve(t(X[,S]) %*% X[, S], t(X[,S]) %*% y)
            beta[S] = beta_1_I
         }
         else{
            break
         }
      }
   }
   return(beta)
}

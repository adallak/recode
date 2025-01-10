#' Compare graphs
#' Part of this function is borrowed from `pcalg` package.
#'
#' @param est  estimated weighted adjacency matrix
#' @param true  true weighted adjacency matrix
#' @param verbose  display the result; default is false
#' @return f1, tpr, fpr, and tdr
#' @export
comparegraph <- function(est, true, verbose = FALSE) {
   est[abs(est) > 0] = 1
   true[abs(true) > 0] = 1
   p = nrow(est)
   #   diag(est) = 0
   #   diag(true) = 0
   diffm = est - true
   nmbtruegaps = (sum(true == 0) - p)/2
   if (nmbtruegaps == 0) {
      fpr = 1
   } else {
      fpr = (sum(diffm >0) / 2) / nmbtruegaps
   }
   diffm2 = true - est
   nmbtrueedge = sum(true == 1) / 2
   if (nmbtrueedge == 0) {
      tpr = 0
   }else {
      tpr = 1 - (sum(diffm2 > 0)/2) / nmbtrueedge
   }
   trueestedge = nmbtrueedge  - sum(diffm2 > 0)/2
   if (nmbtrueedge == 0) {
      if (trueestedge == 0){
         tdr = 1
      }else {
         tdr = 0
      }
   } else{
      tdr = trueestedge / (sum(est == 1)/2)
   }
   fp = sum(diffm > 0)/2
   fn = sum(diffm2 > 0)/2
   tp = sum(true)/2 - fn
   f1 = (2 * tp) / (2 * tp + fp + fn)
   p  = tp / (tp + fp)
   r  = tp / (tp + fn)
   ret = matrix(0, nrow = 1, ncol = 4)
   colnames(ret) = c("F1","TPR","FPR","TDR")
   ret[1,] = round(c(f1, tpr, fpr, tdr),3)
   if (verbose){
      print(ret)
   }
   return(ret)
}


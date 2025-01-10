#' Generates matrix of F-statistics for cardiodata example in Chapter 3
#' Section 2
#'
#'
#' @param data  cardiodata
#' @param top  select the highest `top` number of F-stat; default is 50
#' @param last select the lowest `last` number of F-stat; default is 100
#'
#' @return
#' *`df` generated F-stat
#' @export
generate_cardio_fstat <- function(data, top = 50, last = 100)
{
      cardiodata = data
      top = top
      last = last
      index_0 = which(colnames(cardiodata) == "0")
      index_1 = which(colnames(cardiodata) == "1")

      ## Estimate the number of observation in each class
      n_0 = length(index_0)
      n_1 = length(index_1)

      ## Calculate the dataset mean,mean and variance for each class
      x_bar = rowMeans(cardiodata)
      x_bar0 = rowMeans(cardiodata[, index_0])
      x_bar1 = rowMeans(cardiodata[, index_1])
      var_0 = apply(cardiodata[,index_0], 1, var)
      var_1 = apply(cardiodata[,index_1], 1, var)

      ## Calculate F statistics
      F_stat = (n_0 * (x_bar0 - x_bar)^2 +  n_1 * (x_bar1 - x_bar)^2) / (
            (n_0 - 1) * var_0 + (n_1 - 1) * var_1)

      ## Choose the highest 50 and lowest 100
      ordered_F = order(F_stat, decreasing = TRUE)
      top50 = ordered_F[1:top]
      last100 = ordered_F[(length(ordered_F) - (last - 1)):length(ordered_F)]
      df = t(cardiodata[c(top50, last100),])
return(df)
}

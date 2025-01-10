############################################################################
### Bias example

#' Plots the bias for Lasso, and MCP as described in Chapter 1
#' Section 4.2
#'
#'
#' @param n.list  list of n values
#' @param beta  coefficient matrix
#' @param nsim number of simulations
#' @param seed random seed
#'
#' @return
#' *`bias.lasso` - estimated bias for Lasso
#' *`bias.mcp`     - estimated bias for MCP
#' @export
bias.plot <- function(n.list, beta, nsim = 100, seed = 1234)
{
      p = nrow(beta)
      n.list = n.list
      n.size = length(n.list)
      beta = beta
      ## Matrix to store bias
      bias.lasso = matrix(NA, nrow = n.size, ncol = 4)
      bias.mcp = matrix(NA, nrow = n.size, ncol = 4)
      progress_bar = txtProgressBar(min=0, max=nsim * n.size,
                                    style = 3, char="=")

      nsim = nsim
      n.iter = 0
      for (i in 1:n.size){
            n = n.list[i]
            ## Matrix to store estimated beta
            beta.lasso = matrix(NA, nrow = nsim, ncol = 4)
            beta.mcp = matrix(NA, nrow = nsim, ncol = 4)
            for (sim in 1:nsim){
                  n.iter = n.iter + 1
                  setTxtProgressBar(progress_bar, value = n.iter)
                  set.seed(seed*sim)
                  X = matrix(rnorm(n = n * p), n, p)
                  y = X %*% beta + rnorm(n)
                  ## Implement Lasso
                  cvfit=glmnet::cv.glmnet(x=X, y=y)
                  coef=coef(cvfit,s="lambda.min")
                  beta.lasso[sim,] = coef[2:5]

                  ## Implement MCP
                  cvfit=ncvreg::cv.ncvreg(X=X, y=y, penalty = "MCP")
                  coef=coef(cvfit,s="lambda.min")
                  beta.mcp[sim,] = coef[2:5]
            }
            bias.lasso[i,] = colMeans(beta.lasso) - beta[1:4]
            bias.mcp[i,] = colMeans(beta.mcp) - beta[1:4]
      }

      bias.lasso.df = data.frame(bias.lasso)
      bias.lasso.df["sample"] = n.list

      bias.mcp.df = data.frame(bias.mcp)
      bias.mcp.df["sample"] = n.list

      colors <- c("Unbiased" = "black", "Lasso" = "red", "MCP" = "blue")

      first_coef = ggplot2::ggplot(data=bias.lasso.df, ggplot2::aes(x=sample, y=X1, group=1)) +
         ggplot2::geom_line(size= 1,ggplot2::aes(color = "Lasso"))  +
         ggplot2::xlab("Sample") +
         ggplot2::ylab("Bias") +
         ggplot2::ggtitle(expression(beta[1])) +
         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                        legend.title = ggplot2::element_blank()) +
         ggplot2::geom_line(size= 1,ggplot2::aes(y=0,color = "Unbiased")) +
         ggplot2::geom_line(data=bias.mcp.df, ggplot2::aes(x=sample, y=X1, group=1,color = "MCP")) +
         ggplot2::scale_color_manual(values = colors)

      second_coef = ggplot2::ggplot(data=bias.lasso.df, ggplot2::aes(x=sample, y=X2, group=1)) +
         ggplot2::geom_line(size= 1,ggplot2::aes(color = "Lasso"))  +
         ggplot2::xlab("Sample") + ggplot2::ylab("Bias") +
         ggplot2::ggtitle(expression(beta[2])) +
         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                        legend.title = ggplot2::element_blank()) +
         ggplot2::geom_line(size= 1,ggplot2::aes(y=0,color = "Unbiased")) +
         ggplot2::geom_line(data=bias.mcp.df, ggplot2::aes(x=sample, y=X2, group=1,color = "MCP")) +
         ggplot2::scale_color_manual(values = colors)

      third_coef = ggplot2::ggplot(data=bias.lasso.df, ggplot2::aes(x=sample, y=X3, group=1)) +
         ggplot2::geom_line(size= 1,ggplot2::aes(color = "Lasso"))  +
         ggplot2::xlab("Sample") + ggplot2::ylab("Bias") +
         ggplot2::ggtitle(expression(beta[3])) +
         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                        legend.title = ggplot2::element_blank()) +
         ggplot2::geom_line(size= 1,ggplot2::aes(y=0,color = "Unbiased")) +
         ggplot2::geom_line(data=bias.mcp.df, ggplot2::aes(x=sample, y=X3, group=1,color = "MCP")) +
         ggplot2::scale_color_manual(values = colors)

      fourth_coef = ggplot2::ggplot(data=bias.lasso.df, ggplot2::aes(x=sample, y=X4, group=1)) +
         ggplot2::geom_line(size= 1,ggplot2::aes(color = "Lasso"))  +
         ggplot2::xlab("Sample") + ggplot2::ylab("Bias") +
         ggplot2::ggtitle(expression(beta[4])) +
         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                        legend.title = ggplot2::element_blank()) +
         ggplot2::geom_line(size= 1,ggplot2::aes(y=0,color = "Unbiased")) +
         ggplot2::geom_line(data=bias.mcp.df, ggplot2::aes(x=sample, y=X4, group=1,color = "MCP")) +
         ggplot2::scale_color_manual(values = colors)

      ggpubr::ggarrange(first_coef, second_coef,
                third_coef,fourth_coef,  common.legend = TRUE,
                legend = "bottom", ncol = 2, nrow = 2)

     # return(list("bias.lasso" = bias.lasso, "bias.mcp" = bias.mcp))
}

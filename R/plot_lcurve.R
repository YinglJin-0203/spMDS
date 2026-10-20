#' Plot the L-curve (log-log) with lambda labels
#'
#' @param lcurve_df output of compute_lcurve()
#' @param label_points if TRUE, annotate each point with its lambda value
plot_lcurve <- function(lcurve_df, label_points = TRUE) {
  plot(lcurve_df$penalty, lcurve_df$stress,
       type = "b", pch = 19, col = "steelblue",
       xlab = "regularization",
       ylab = "stress",
       main = "L-curve for spMDS")
  if (label_points) {
    text(lcurve_df$penalty, lcurve_df$stress,
         labels = signif(lcurve_df$lambda, 3),
         pos = 3, cex = 0.7, col = "gray30")
  }
}

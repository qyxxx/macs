print.forecast.ctree <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  predictor <- x$predictor
  names(predictor) <- 1:length(predictor)
  cat("\n")
  print.default(format(predictor), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

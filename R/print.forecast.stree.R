print.forecast.stree <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  predictor <- x$predictor
  names(predictor) <- 1:length(predictor)
  cat("\n")
  print.default(format(predictor, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

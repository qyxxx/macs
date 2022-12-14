print.forecast.masal <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\n")
  print.default(format(x, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

print.masal <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\ncall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nfitted:\n")
  print.default(format(x$fitted), print.gap = 2L, quote = FALSE)
  cat("\n")

  invisible(x)
}

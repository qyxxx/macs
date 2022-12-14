print.ctree <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\ncall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nnodeclass:\n")
  print.default(format(x$nodeclass), print.gap = 2L, quote = FALSE)
  cat("\nnnd:\n", x$nnd, "\n", sep = "")
  cat("\n")

  invisible(x)
}

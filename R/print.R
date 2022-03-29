# Print contents of \code{rfvimptest} object.
#' @author Alexander Hapfelmeier, Roman Hornung
#' @export
print.rfvimptest <- function(x, ...) {
  cat("\n")
  cat("\t Result of (sequential) permutation tests for statistical significance of predictors in random forests\n")
  cat("\n")
  cat("Type of test:                    ", x$testtype, "\n")
  cat("Number of tested variables:      ", length(x$varimp), "\n")
  cat("Number (proportion) significant: ", paste0(sum(x$testres=="accept H1"), " (", mean(x$testres=="accept H1"), ")"), "\n")
  cat("\n")
  cat("Number of trees:                 ", x$ntree, "\n")
  cat("Maximum number of permutations:  ", x$Mmax, "\n")
  cat("Computation time:                ", x$comptime, "\n\n")
}

# Print contents of \code{allinone} object.
#' @author Alexander Hapfelmeier, Roman Hornung
#' @export
print.allinone <- function(x, ...) {
  cat("\n")
  cat("\t Comparison of all applicable (sequential) permutation testing approaches\n")
  cat("\n")
  cat("Compared approaches: ", paste(names(x$testres)[-1], collapse=", "), "\n")
  cat("\n")
  cat("Number of tested variables:      ", length(x$varimp), "\n")
  cat("\n")
  cat("Number of trees:                 ", x$ntree, "\n")
  cat("Maximum number of permutations:  ", x$Mmax, "\n")
  cat("Computation time:                ", x$comptime, "\n\n")
}

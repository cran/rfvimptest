##' Threshold values of a SPRT and SAPT, and calculation of sequential Monte Carlo p-values
##'
##' Function to produce threshold values for the number of cumulative successes a sequential permutation test of \code{type="SRPT"} or \code{type="SAPT"} needs to show in order to "keep H0" or "accept H1" early in the sequence of tests.
##'
##' An example of how to calculate sequential Monte Carlo p-values is also given.
##'
##' The function and examples are intended to support the use of SPRT, SAPT and sequential p-values for any type of modelling approach, even outside the application in \code{rfvimptest::rfvimptest()}.
##'
##' @param p0 The value of the p-value in the null hypothesis (H0: p = p0) of SPRT and SAPT. Default is 0.06.
##' @param p1 The value of the p-value in the alternative hypothesis (H1: p = p1) of SPRT and SAPT. Default is 0.04.
##' @param alpha The significance level of SPRT when p = p0. Also known as type I error. Default is 0.05.
##' @param beta One minus the power of SPRT when p = p1. Also known as type II error. Default is 0.2.
##' @param A The quantity A in the formula of SAPT. Default is 0.1 for a type I error of 0.05. Usually not changed by the user.
##' @param B The quantity B in the formula of SAPT. Default is 10 (1/A) for a type I error of 0.05. Usually not changed by the user.
##' @param Mmax Maximum number of permutations used in each permutation test. Default is 500.
##' @param type Type of the sequential test method. The choices are "SPRT" and "SAPT".
##' @return List of two vectors, each containing \code{Mmax} threshold values for the number of cumulative successes a sequential permutation test of \code{type="SRPT"} or \code{type="SAPT"} needs to show in order to "keep H0" or "accept H1" early in the sequence of tests.
##' @examples
##' \donttest{
##'
##' ## Load package:
##' library("rfvimptest")
##'
##' ## Set seed to obtain reproducible results:
##' set.seed(1234)
##'
##'
##' # Load example data
##' data(hearth2)
##'
##'
##' # Calculation of a SAPT
##'
##' # Create critical values for a SAPT
##' myvals <- threshold_values(p0 = 0.06, p1 = 0.04, A = 0.1, B = 10, Mmax = 200, type = "SAPT")
##' myvals
##'
##' # Fit a model to the original data
##' mod <- glm(I(Class == 2) ~ ., data = hearth2)
##' summary(mod)
##'
##' # Derive a statistic of interest from the model. Here, the negative AIC is used as a statistic
##' # to be maximised.
##' stat <- -1 * mod$aic
##'
##' # Perform the permutations, extract the statistics and stop early according to the critical values
##' myresult <- sapply(setdiff(names(hearth2), "Class"), function(j) {
##'   permdat <- hearth2
##'   permstats <- c()
##'   for (i in 1:length(myvals[[1]])) {
##'     permdat[, j] <- sample(permdat[, j])
##'     permstats <- c(permstats, -1 * glm(I(Class == 2) ~ ., data = permdat)$aic)
##'     if (sum(permstats >= stat) >= myvals[[1]][i] | sum(permstats >= stat) <= myvals[[2]][i]) break
##'   }
##'   permstats
##' })
##'
##' # Statistics obtained after permutation of the data
##' myresult
##'
##' # Number of permutations performed until (early) stopping
##' sapply(myresult, length)
##'
##' # Derive the SAPT decisions from the statistics
##' sapply(myresult, function(x) {
##'   if (sum(x >= stat) >= myvals[[1]][length(x)]) return("keep H0")
##'   else if (sum(x >= stat) <= myvals[[2]][length(x)]) return("accept H1")
##'   else ifelse(sum(x >= stat) / length(myvals[[1]]) <= 0.05, "accept H1", "keep H0")
##' })
##'
##'
##' # Calculation of sequential Monte Carlo p-values
##'
##' h <- 10 # Parameter h of the sequential Monte Carlo p-value calculation.
##'
##' # Perform the permutations, extract the statistics and stop early when h successes are reached
##' mypvals <- sapply(setdiff(names(hearth2), "Class"), function(j) {
##'   permdat <- hearth2
##'   permstats <- c()
##'   Mmax <- length(myvals[[1]])
##'   for (i in 1:Mmax) {
##'     permdat[, j] <- sample(permdat[, j])
##'     permstats <- c(permstats, -1 * glm(I(Class == 2) ~ ., data = permdat)$aic)
##'     d <- sum(permstats >= stat)
##'     if (d == h) {pval <- d/length(permstats); m <- i; break}
##'     else if (i == Mmax) {pval <- (d + 1)/(Mmax + 1); m <- i}
##'   }
##'   c(pval, m)
##' })
##'
##' # p-values and number of permutations performed until (early) stopping
##' rownames(mypvals) <- c("p-value", "m")
##' t(mypvals)
##'
##' }
##'
##' @author Alexander Hapfelmeier, Roman Hornung
##' @references
##' \itemize{
##'   \item Wald, A. (1945). Sequential tests of statistical hypotheses. Ann. Math. Stat. 16, 117-186, <\doi{10.1214/aoms/1177731118}>.
##'   \item Lock, R.H. (1991). A sequential approximation to a permutation test. Commun. Stat., Simul. Comput. 20, 341-363, <\doi{10.1080/03610919108812956}>.
##'   \item Besag, J., Clifford, P. (1991). Sequential Monte Carlo p-values. Biometrika 78, 301-304, <\doi{10.1093/biomet/78.2.301}>.
##'   \item Hapfelmeier, A., Hornung, R., Haller, B. (2023). Efficient permutation testing of variable importance measures by the example of random forests. Comput Stat Data Anal, 181:107689, <\doi{10.1016/j.csda.2022.107689}>.
##'   }
##' @encoding UTF-8
##' @export
threshold_values <- function(p0 = 0.06, p1 = 0.04, alpha = 0.05, beta = 0.2, A = 0.1, B = 10, Mmax = 500, type = c("SPRT", "SAPT")[1]) {
  if (type == "SPRT") {
    A <- beta / (1 - alpha)
    B <- (1 - beta) / alpha
  }

  logA <- log(A)
  logB <- log(B)
  help1 <- log((1 - p0) / (1 - p1))
  help2 <- log((p1 * (1 - p0)) / (p0 * (1 - p1)))

  list((logA + 1:Mmax * help1) / help2, (logB + 1:Mmax * help1) / help2)
}

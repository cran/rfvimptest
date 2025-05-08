##' Apply all available (sequential) permutation testing approaches of variable importance measures with one function call
##'
##' This is a helper function, which allows to perform all (sequential) permutation testing approaches of variable importance measures described in \code{\link{rfvimptest}}
##' with a single function call. This may be useful for comparing the results obtained using the different approaches.
##' Importantly, this function is computationally efficient by re-using the permuted variable importance values obtained
##' for the conventional permutation test (that performs all \code{Mmax} permutations) for the other approaches. For details
##' on the different approaches see \code{\link{rfvimptest}}.
##'
##' @param data A \code{data.frame} containing the variables in the model.
##' @param yname Name of outcome variable.
##' @param Mmax Maximum number of permutations used in each permutation test. Default is 500.
##' @param varnames Optional. Names of the variables for which testing should be performed. By default all variables in \code{data} with the exception of the outcome variable are used.
##' @param p0 The value of the p-value in the null hypothesis (H0: p = p0) of SPRT and SAPT. Default is 0.06.
##' @param p1 The value of the p-value in the alternative hypothesis (H1: p = p1) of SPRT and SAPT. Default is 0.04.
##' @param alpha The significance level of SPRT when p = p0. Also known as type I error. Default is 0.05.
##' @param beta One minus the power of SPRT when p = p1. Also known as type II error. Default is 0.2.
##' @param A The quantity A in the formula of SAPT. Default is 0.1 for a type I error of 0.05. Usually not changed by the user.
##' @param B The quantity B in the formula of SAPT. Default is 10 (1/A) for a type I error of 0.05. Usually not changed by the user.
##' @param h The quantity h in the formula for the sequential Monte Carlo p-value. The default value for h is 8. Larger values lead to more precise p-value estimates,
##' but are computationally more expensive.
##' @param nperm The numbers of permutations of the out-of-bag observations over which the results are averaged, when calculating the variable importance measure values. Default is 1. Larger values than 1 can only be considered when \code{condinf=TRUE}, that is, when using random forests
##' with conditional inference trees (Hothorn et al., 2006) as base learners.
##' @param ntree Number of trees per forest. Default is 500.
##' @param progressbar Output the current progress of the calculations for each variable to the console? Default is TRUE.
##' @param condinf Set this value to \code{TRUE} if random forests using conditional inference trees (Hothorn et al., 2006) should
##' be used and to \code{FALSE} if classical random forests using CART trees should be used. Default is \code{FALSE}.
##' @param ... Further arguments passed to \code{ranger::ranger} (if \code{condinf=FALSE}) or \cr \code{party::cforest_unbiased()} (if \code{condinf=TRUE}).
##' @return Object of class \code{allinone} with elements
##'   \item{\code{varimp}}{Variable importance for each considered independent variable.}
##'   \item{\code{testres}}{The results ("keep H0" vs. "accept H1") of the tests for each considered independent variable.}
##'   \item{\code{pvalues}}{The p-values of the tests for each considered independent variable. Note that p-values are only obtained for the
##'   method types "pval" and "complete".}
##'   \item{\code{stoppedearly}}{For each independent variable, whether the calculations stopped early ("yes") or the maximum of \code{Mmax} permutations was reached ("no").}
##'   \item{\code{perms}}{The number of permutations performed for each independent variable.}
##'   \item{\code{Mmax}}{Maximum number of permutations used in each permutation test.}
##'   \item{\code{ntree}}{Number of trees per forest.}
##'   \item{\code{comptime}}{The time the computations needed.}
##' @examples
##' \donttest{
##'
##' # Load package:
##' library("rfvimptest")
##'
##' # Set seed to obtain reproducible results:
##' set.seed(1234)
##'
##' # Load example data:
##' data(hearth2)
##'
##' # NOTE: For illustration purposes very small numbers of maximum
##' # permutations are considered in the below examples.
##' # This number would be much too small for actual applications.
##' # The default number is Max=500.
##'
##' # When using condinf=FALSE (default) the results for the two-sample
##' # permutation tests are not obtained:
##' (ptest <- allinone(data=hearth2, yname="Class",  Mmax=20))
##'
##' # Variable importance values with p-values from the Monte Carlo p-value
##' # and the complete approach:
##' ptest$varimp
##' ptest$pvalues$pval
##' ptest$pvalues$complete
##'
##'
##' # When setting condinf=TRUE the results are obtained for all approaches,
##' # that is, including those for the two-sample permutation tests
##' # (in this illustration very small number of trees ntree=30 are used,
##' # in practice much larger numbers should be used; the default is ntree=500):
##' (ptest_ci <- allinone(data=hearth2, yname="Class", condinf=TRUE, ntree=30, Mmax=10))
##' ptest_ci$testres
##'
##' }
##'
##' @author Alexander Hapfelmeier, Roman Hornung
##' @references
##' \itemize{
##'   \item Breiman, L. (2001). Random forests. Mach Learn, 45:5-32, <\doi{10.1023/A:1010933404324}>.
##'   \item Coleman, T., Peng, W., Mentch, L. (2019). Scalable and efficient hypothesis testing with random forests. arXiv preprint arXiv:1904.07830, <\doi{10.48550/arXiv.1904.07830}>.
##'   \item Hapfelmeier, A., Hornung, R., Haller, B. (2023). Efficient permutation testing of variable importance measures by the example of random forests. Comput Stat Data Anal, 181:107689, <\doi{10.1016/j.csda.2022.107689}>.
##'   \item Hapfelmeier, A., Ulm, K. (2013). A new variable selection approach using Random Forests. CSDA 60:50–69, <\doi{10.1016/j.csda.2012.09.020}>.
##'   \item Hapfelmeier, A., Hothorn, T., Ulm, K., Strobl, C. (2014). A new variable importance measure for random forests with missing data. Stat Comput 24:21–34, <\doi{10.1007/s11222-012-9349-1}>.
##'   \item Hothorn, T., Hornik, K., Zeileis, A. (2006). Unbiased Recursive Partitioning: A Conditional Inference Framework. J Comput Graph Stat 15(3):651–674, <\doi{10.1198/106186006X133933}>.
##'   \item Wright, M. N., Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17, <\doi{10.18637/jss.v077.i01}>.
##'   }
##' @seealso \code{\link{rfvimptest}}
##' @encoding UTF-8
##' @importFrom stats as.formula
##' @importFrom utils setTxtProgressBar txtProgressBar
##' @export
allinone <- function(data, yname, Mmax = 500, varnames = NULL, p0 = 0.06, p1 = 0.04, alpha = 0.05, beta = 0.2, A = 0.1, B = 10, h = 8, nperm = 1,
                     ntree = 500,
                     progressbar = TRUE, condinf=FALSE, ...) {
  starttime <- Sys.time()

  if(any(is.na(data))) {
    missvariables <- paste(names(data)[apply(data, 2, function(x) any(is.na(x)))], collapse = ", ")
    stop(paste0("Missing data in columns: ", missvariables, ". Please provide complete data or consider setting condinf=TRUE."))
  }

  if(!condinf)
    warning("Results for the 'twosample' approach were not calculated because condinf=FALSE.")

  if(!condinf & nperm > 1)
    stop("Values of 'nperm' different than 1 are only useable if condinf=TRUE.")

  if (progressbar) pb <- txtProgressBar(min = 1, max = Mmax, initial = 1, width = 10, style = 3, char = "|")

  # Define the parameters A and B for SPRT and SAPT
  A_SPRT <- beta / (1 - alpha)
  B_SPRT <- (1 - beta) / alpha

  logA_SPRT <- log(A_SPRT)
  logB_SPRT <- log(B_SPRT)

  logA_SAPT <- log(A)
  logB_SAPT <- log(B)

  # Define the stopping criteria
  help1 <- log((1 - p0) / (1 - p1))
  help2 <- log((p1 * (1 - p0)) / (p0 * (1 - p1)))

  stop_crits <- list(SPRT = list((logA_SPRT + 1:Mmax * help1) / help2, (logB_SPRT + 1:Mmax * help1) / help2),
                     SAPT = list((logA_SAPT + 1:Mmax * help1) / help2, (logB_SAPT + 1:Mmax * help1) / help2),
                     certain = list(rep(alpha*Mmax, times = Mmax), Mmax*alpha - Mmax + 1:Mmax))

  if (is.null(varnames))
    varnames <- names(data)[names(data)!=yname]

  if (!condinf) {
    rfmod <- ranger::ranger(data = data, dependent.variable.name=yname, num.tree=ntree, importance = "permutation")
    vimp_orig <- list()
    vimp_orig$values <- rfmod$variable.importance[varnames]
  }
  else {
    rfmod <- party::cforest(as.formula(paste(yname, " ~ .", sep="")), data = data, controls = party::cforest_unbiased(ntree=ntree, ...))
    vimp_orig <- permimp::permimp(rfmod, whichxnames = varnames, nperm = nperm, asParty = TRUE, progressBar = FALSE)
  }



  # The computations for each variable

  methodnames <- c("SPRT", "SAPT", "pval", "certain", "complete", "SPRT_TS", "SAPT_TS", "pval_TS", "certain_TS", "complete_TS")
  methodnamespval <- c("pval", "complete", "pval_TS", "complete_TS")

  testresult <- lapply(varnames, function(v) {
    # copy data for permutation
    permdata <- data
    permvimps <- rep(NA, Mmax)

    # do the permutations of the gereral permutation test
    for (m in 1:Mmax) {
      if (progressbar) {setTxtProgressBar(pb, m)}
      if (m == 1 & progressbar) cat(" of variable", v)
      permdata[, v] <- sample(permdata[, v])
      if (!condinf) {
        permmod <- ranger::ranger(data = permdata, dependent.variable.name=yname, num.tree=ntree, importance="permutation")
        permvimps[m] <- permmod$variable.importance[v]
      }
      else {
        permmod <- party::cforest(as.formula(paste(yname, " ~ .", sep="")), data = permdata, controls = party::cforest_unbiased(ntree=ntree, ...))
        permmodvimps <- permimp::permimp(permmod, whichxnames = v, nperm = nperm, asParty = TRUE, progressBar = FALSE)
        permvimps[m] <- permmodvimps$values
      }
    }

    # compute d along the sequence 1:Mmax and the p-value of COMPLETE
    d <- cumsum(permvimps >= vimp_orig$values[v])
    pval <- rev(d)[1] / Mmax

    # apply the stopping criteria for each method to d
    m1_SPRT <- which(d >= stop_crits[[1]][[1]])[1]
    m2_SPRT <- which(d <= stop_crits[[1]][[2]])[1]

    # if a method does not stop early m is set to Mmax+1 for ease of further processing
    if (is.na(m1_SPRT)) m1_SPRT <- Mmax + 1
    if (is.na(m2_SPRT)) m2_SPRT <- Mmax + 1

    m1_SAPT <- which(d >= stop_crits[[2]][[1]])[1]
    m2_SAPT <- which(d <= stop_crits[[2]][[2]])[1]

    if (is.na(m1_SAPT)) m1_SAPT <- Mmax + 1
    if (is.na(m2_SAPT)) m2_SAPT <- Mmax + 1

    m_PVAL <- which(d == h)[1]

    if (is.na(m_PVAL)) m_PVAL <- Mmax + 1

    m1_CERTAIN <- which(d > stop_crits[[3]][[1]])[1]
    m2_CERTAIN <- which(d <= stop_crits[[3]][[2]])[1]

    if (is.na(m1_CERTAIN)) m1_CERTAIN <- Mmax + 1
    if (is.na(m2_CERTAIN)) m2_CERTAIN <- Mmax + 1

    if (condinf) {

      # Repeat all steps for the two-sample test
      permvimps_TS <- rep(NA, Mmax)

      for (m in 1:Mmax) {
        permvimps_TS[m] <- mean(c(vimp_orig$perTree[, v], permmodvimps$perTree[, v])[sample(1:(2*ntree), ntree)])
      }

      d_TS <- cumsum(permvimps_TS >= vimp_orig$values[v])
      pval_TS <- rev(d_TS)[1] / Mmax

      m1_SPRT_TS <- which(d_TS >= stop_crits[[1]][[1]])[1]
      m2_SPRT_TS <- which(d_TS <= stop_crits[[1]][[2]])[1]

      if (is.na(m1_SPRT_TS)) m1_SPRT_TS <- Mmax + 1
      if (is.na(m2_SPRT_TS)) m2_SPRT_TS <- Mmax + 1

      m1_SAPT_TS <- which(d_TS >= stop_crits[[2]][[1]])[1]
      m2_SAPT_TS <- which(d_TS <= stop_crits[[2]][[2]])[1]

      if (is.na(m1_SAPT_TS)) m1_SAPT_TS <- Mmax + 1
      if (is.na(m2_SAPT_TS)) m2_SAPT_TS <- Mmax + 1

      m_PVAL_TS <- which(d_TS == h)[1]

      if (is.na(m_PVAL_TS)) m_PVAL_TS <- Mmax + 1

      m1_CERTAIN_TS <- which(d_TS > stop_crits[[3]][[1]])[1]
      m2_CERTAIN_TS <- which(d_TS <= stop_crits[[3]][[2]])[1]

      if (is.na(m1_CERTAIN_TS)) m1_CERTAIN_TS <- Mmax + 1
      if (is.na(m2_CERTAIN_TS)) m2_CERTAIN_TS <- Mmax + 1

      testres <- stoppedearly <- permperf <- rep(NA, 10)
      pvalue <- rep(NA, 4)

    }
    else {
      testres <- stoppedearly <- permperf <- rep(NA, 5)
      pvalue <- rep(NA, 2)
    }

    names(testres) <- names(stoppedearly) <- names(permperf) <- methodnames[1:length(testres)]
    names(pvalue) <- methodnamespval[1:length(pvalue)]


    if (any(c(m1_SPRT, m2_SPRT) < Mmax)) {
      testres[1] <- c("keep H0", "accept H1")[which.min(c(m1_SPRT, m2_SPRT))]
      permperf[1] <- c(m1_SPRT, m2_SPRT)[which.min(c(m1_SPRT, m2_SPRT))]
      stoppedearly[1] <- "yes"
    }
    else {
      testres[1] <- ifelse(pval > 0.05, "keep H0", "accept H1")
      permperf[1] <- Mmax
      stoppedearly[1] <- "no"
    }


    if (any(c(m1_SAPT, m2_SAPT) < Mmax)) {
      testres[2] <- c("keep H0", "accept H1")[which.min(c(m1_SAPT, m2_SAPT))]
      permperf[2] <- c(m1_SAPT, m2_SAPT)[which.min(c(m1_SAPT, m2_SAPT))]
      stoppedearly[2] <- "yes"
    }
    else {
      testres[2] <- ifelse(pval > 0.05, "keep H0", "accept H1")
      permperf[2] <- Mmax
      stoppedearly[2] <- "no"
    }


    if (m_PVAL <= Mmax) {
      pvalue[1] <- h / m_PVAL
      testres[3] <- ifelse(pvalue[1] > 0.05, "keep H0", "accept H1")
      permperf[3] <- m_PVAL
      stoppedearly[3] <- "yes"
    }
    else {
      pvalue[1] <- (rev(d)[1] + 1) / (Mmax + 1)
      testres[3] <- ifelse(pvalue[1] > 0.05, "keep H0", "accept H1")
      permperf[3] <- Mmax
      stoppedearly[3] <- "no"
    }


    if (any(c(m1_CERTAIN, m2_CERTAIN) < Mmax)) {
      testres[4] <- c("keep H0", "accept H1")[which.min(c(m1_CERTAIN, m2_CERTAIN))]
      permperf[4] <- c(m1_CERTAIN, m2_CERTAIN)[which.min(c(m1_CERTAIN, m2_CERTAIN))]
      stoppedearly[4] <- "yes"
    }
    else {
      testres[4] <- ifelse(pval > 0.05, "keep H0", "accept H1")
      permperf[4] <- Mmax
      stoppedearly[4] <- "no"
    }


    pvalue[2] <- pval
    testres[5] <- ifelse(pvalue[2] > 0.05, "keep H0", "accept H1")
    permperf[5] <- Mmax
    stoppedearly[5] <- "no"


    if (condinf) {


      if (any(c(m1_SPRT_TS, m2_SPRT_TS) < Mmax)) {
        testres[6] <- c("keep H0", "accept H1")[which.min(c(m1_SPRT_TS, m2_SPRT_TS))]
        permperf[6] <- c(m1_SPRT_TS, m2_SPRT_TS)[which.min(c(m1_SPRT_TS, m2_SPRT_TS))]
        stoppedearly[6] <- "yes"
      }
      else {
        testres[6] <- ifelse(pval_TS > 0.05, "keep H0", "accept H1")
        permperf[6] <- Mmax
        stoppedearly[6] <- "no"
      }


      if (any(c(m1_SAPT_TS, m2_SAPT_TS) < Mmax)) {
        testres[7] <- c("keep H0", "accept H1")[which.min(c(m1_SAPT_TS, m2_SAPT_TS))]
        permperf[7] <- c(m1_SAPT_TS, m2_SAPT_TS)[which.min(c(m1_SAPT_TS, m2_SAPT_TS))]
        stoppedearly[7] <- "yes"
      }
      else {
        testres[7] <- ifelse(pval_TS > 0.05, "keep H0", "accept H1")
        permperf[7] <- Mmax
        stoppedearly[7] <- "no"
      }


      if (m_PVAL_TS <= Mmax) {
        pvalue[3] <- h / m_PVAL_TS
        testres[8] <- ifelse(pvalue[3] > 0.05, "keep H0", "accept H1")
        permperf[8] <- m_PVAL_TS
        stoppedearly[8] <- "yes"
      }
      else {
        pvalue[3] <- (rev(d_TS)[1] + 1) / (Mmax + 1)
        testres[8] <- ifelse(pvalue[3] > 0.05, "keep H0", "accept H1")
        permperf[8] <- Mmax
        stoppedearly[8] <- "no"
      }


      if (any(c(m1_CERTAIN_TS, m2_CERTAIN_TS) < Mmax)) {
        testres[9] <- c("keep H0", "accept H1")[which.min(c(m1_CERTAIN_TS, m2_CERTAIN_TS))]
        permperf[9] <- c(m1_CERTAIN_TS, m2_CERTAIN_TS)[which.min(c(m1_CERTAIN_TS, m2_CERTAIN_TS))]
        stoppedearly[9] <- "yes"
      }
      else {
        testres[9] <- ifelse(pval_TS > 0.05, "keep H0", "accept H1")
        permperf[9] <- Mmax
        stoppedearly[9] <- "no"
      }


      pvalue[4] <- pval_TS
      testres[10] <- ifelse(pvalue[4] > 0.05, "keep H0", "accept H1")
      permperf[10] <- Mmax
      stoppedearly[10] <- "no"

    }

    if (progressbar) cat(" - finished", "\n")

    return(list(testres=testres, pvalue=pvalue, stoppedearly=stoppedearly, permperf=permperf))

  })

  stoptime <- Sys.time()
  time_elapsed <- paste0(round(as.numeric(difftime(stoptime, starttime, units = "secs")), 1), " seconds")

  testres <- as.data.frame(t(sapply(testresult, function(x) x$testres)))
  pvalues <- as.data.frame(t(sapply(testresult, function(x) x$pvalue)))
  stoppedearly <- as.data.frame(t(sapply(testresult, function(x) x$stoppedearly)))
  perms <- as.data.frame(t(sapply(testresult, function(x) x$permperf)))

  testres$variable <- pvalues$variable <- stoppedearly$variable <- perms$variable <- names(vimp_orig$values)
  testres <- testres[,c(ncol(testres), 1:(ncol(testres)-1))]
  pvalues <- pvalues[,c(ncol(pvalues), 1:(ncol(pvalues)-1))]
  stoppedearly <- stoppedearly[,c(ncol(stoppedearly), 1:(ncol(stoppedearly)-1))]
  perms <- perms[,c(ncol(perms), 1:(ncol(perms)-1))]

  names(testres)[1] <- names(pvalues)[1] <- names(stoppedearly)[1] <- names(perms)[1] <- "variable"

  result <- list(varimp=vimp_orig$values, testres = testres, pvalues = pvalues,
                 stoppedearly = stoppedearly, perms = perms, Mmax=Mmax, ntree=ntree, comptime = time_elapsed)

  class(result) <- "allinone"

  return(result)

}

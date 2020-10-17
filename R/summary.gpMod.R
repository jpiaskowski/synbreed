#' Summary for class gpMod
#'
#' S3 \code{summary} method for objects of class \code{gpMod}
#'
#'
#' @aliases summary.gpMod  summary.gpModList  print.summary.gpMod print.summary.gpModList
#' @param object object of class \code{gpMod}
#' @seealso \code{\link{gpMod}}
#' @rdname gpMod
#' @examples
#'
#' \dontrun{
#' library(synbreedData)
#' data(maize)
#' maizeC <- codeGeno(maize)
#' # marker-based (realized) relationship matrix
#' U <- kin(maizeC, ret = "realized") / 2
#'
#' # BLUP model
#' mod <- gpMod(maizeC, model = "BLUP", kin = U)
#' summary(mod)
#' }
#'
#' @export
summary.gpMod <- function(object, ...) {
  ans <- list()
  ans$model <- object$model
  if (object$model %in% c("BLUP")) ans$summaryFit <- summary(object$fit)
  if (object$model == "BL") ans$summaryFit <- list(mu = object$fit$mu, varE = object$fit$varE, varU = object$fit$varU, lambda = object$fit$lambda, nIter = object$fit$nIter, burnIn = object$fit$burnIn, thin = object$fit$thin)
  if (object$model == "BRR") ans$summaryFit <- list(mu = object$fit$mu, varE = object$fit$varE, varBr = object$fit$varBr, varU = object$fit$varU, nIter = object$fit$nIter, burnIn = object$fit$burnIn, thin = object$fit$thin)
  ans$n <- sum(!is.na(object$y))
  ans$sumNA <- sum(is.na(object$y))
  ans$summaryG <- summary(as.numeric(object$g))
  if (is.null(object$prediction)) ans$summaryP <- NULL else ans$summaryP <- summary(as.numeric(object$prediction))
  class(ans) <- "summary.gpMod"
  ans
}

summary.gpModList <- function(object, ...) {
  ret <- list()
  for (i in 1:length(object)) {
    ret[[i]] <- summary(object[[i]])
  }
  class(ret) <- "summary.gpModList"
  names(ret) <- names(object)
  ret
}

print.summary.gpMod <- function(x, ...) {
  cat("Object of class 'gpMod' \n")
  cat("Model used:", x$model, "\n")
  cat("Nr. observations ", x$n, " \n", sep = "")
  cat("Genetic performances: \n")
  cat("  Min.    1st Qu. Median  Mean    3rd Qu. Max    \n")
  cat(format(x$summaryG, width = 7, trim = TRUE), "\n", sep = " ")
  if (!is.null(x$summaryP)) {
    cat("\nGenetic performances of predicted individuals: \n", sep = " ")
    cat("  Min.    1st Qu. Median  Mean    3rd Qu. Max    \n")
    cat(format(x$summaryP, width = 7, trim = TRUE), "\n", sep = " ")
  }
  cat("--\n")
  cat("Model fit \n")
  if (x$model %in% c("BLUP")) {
    cat(print(x$summaryFit), "\n")
  } else {
    cat("MCMC options: nIter = ", x$summaryFit$nIter, ", burnIn = ", x$summaryFit$burnIn, ", thin = ", x$summaryFit$thin, "\n", sep = "")
    cat("             Posterior mean \n")
    cat("(Intercept) ", x$summaryFit$mu, "\n")
    cat("VarE        ", x$summaryFit$varE, "\n")
    if (!is.null(x$summaryFit$varU)) cat("VarU        ", x$summaryFit$varU, "\n")
    if (x$model == "BL") {
      cat("lambda      ", x$summaryFit$lambda, "\n")
    }
    if (x$model == "BRR") {
      cat("varBr       ", x$summaryFit$varBr, "\n")
    }
  }
}

print.summary.gpModList <- function(x, ...) {
  for (i in 1:length(x)) {
    cat(paste("\n\tTrait ", names(x)[i], "\n\n\n"))
    print(x[[i]])
    if (i != length(x)) cat("-------------------------------------------------\n\n")
  }
}

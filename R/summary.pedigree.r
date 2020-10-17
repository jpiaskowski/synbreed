#' Summary of pedigree information
#'
#' \code{Summary} method for class ''pedigree''
#'
#'
#' @aliases summary.pedigree print.summary.pedigree
#' @param object object of class ''pedigree''
#' @param ...  not used
#' @author Valentin Wimmer
#' @examples
#'
#' # plant pedigree
#' ped <- simul.pedigree(gener = 4, 7)
#' summary(ped)
#'
#' # animal pedigree
#' ped <- simul.pedigree(gener = 4, 7, animals = TRUE)
#' summary(ped)
#' @export
summary.pedigree <- function(object, ...) {
  if (any(class(object) == "gpData")) {
    ped <- object$pedigree
  } else {
    ped <- object
  }

  n <- nrow(ped)
  ans <- list(nID = n, nPar1 = length(unique(ped$Par1[ped$Par1 != 0])), nPar2 = length(unique(ped$Par2[ped$Par2 != 0])), nGener = length(unique(ped$gener)), nUnknownParents = sum(ped$Par1 == 0) + sum(ped$Par2 == 0))
  if (!is.null(ped$sex)) ans$sex <- c(males = sum(ped$sex), females = sum(1 - ped$sex))
  class(ans) <- "summary.pedigree"
  ans
}

# print method for summary.pedigree
print.summary.pedigree <- function(x, ...) {
  cat("Number of \n")
  cat("\t individuals ", x$nID, "\n")
  if (!is.null(x$sex)) {
    cat("\t males : ", x$sex[1], ", females : ", x$sex[2], "\n")
    cat("\t Par 1 (sire) ", x$nPar1, "\n")
    cat("\t Par 2 (dam)  ", x$nPar2, "\n")
  }
  else {
    cat("\t Par 1       ", x$nPar1, "\n")
    cat("\t Par 2       ", x$nPar2, "\n")
  }
  cat("\t generations ", x$nGener, "\n")
}

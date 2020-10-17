#' Prediction for genomic prediction models.
#'
#' S3 \code{predict} method for objects of class \code{gpMod}. A genomic
#' prediction model is used to predict the genetic performance for e.g.
#' unphenotyped individuals using an object of class \code{gpMod} estimated by
#' a training set.
#'
#' For models, \code{model="RR"} and \code{"BL"}, the prediction for the
#' unphenotyped individuals is given by \deqn{\hat{g}'=\hat{\mu} +
#' W'\hat{m}}{ghat=muhat + W mhat} with the estimates taken from the
#' \code{gpMod} object. For the prediction using \code{model="BLUP"}, the full
#' relationship matrix including individuals of the training set and the
#' prediction set must be specified in the \code{gpMod}. This model is used to
#' predict the unphenotyped individuals of the prediction set by solving the
#' corresponding mixed model equations using the variance components of the fit
#' in \code{gpMod}.
#'
#' @aliases predict.gpMod
#'
#' @param object object of class \code{gpMod} which is the model used for the
#' prediction. If the model includes a \code{relationshipMatrix}, this must
#' include both the individuals in the training data used for fitting
#' \code{gpMod} and those which sould be predicted in \code{newdata} (see
#' example below).
#' @param newdata for \code{model="BL"} and \code{"BRR"} an object of class
#' \code{gpData} with the marker data of the unphenotyped individuals. For
#' \code{model="BLUP"} a \code{character} vector with the names of the
#' individuals to predict. If \code{newdata=NULL}, the genetic performances of
#' the individuals for the training set are returned.
#' @param \dots not used
#'
#' @return a named vector with the predicted genetic values for all individuals
#' in \code{newdata}.
#' @author Valentin Wimmer
#' @seealso \code{\link{gpMod}}
#' @references Henderson C (1977) Best linear unbiased prediction of breeding
#' values not in the model for records. Journal of Dairy Science 60:783-787
#'
#' Henderson CR (1984). Applications of linear models in animal breeding.
#' University of Guelph.
#' @examples
#'
#' # Example from Henderson (1977)
#' dat <- data.frame(y = c(132, 147, 156, 172), time = c(1, 2, 1, 2), row.names = c("ID1", "ID2", "ID3", "ID4"))
#' ped <- create.pedigree(
#'   ID = c("ID6", "ID5", "ID1", "ID2", "ID3", "ID4"),
#'   Par1 = c(0, 0, "ID5", "ID5", "ID1", "ID6"),
#'   Par2 = c(0, 0, 0, 0, "ID6", "ID2")
#' )
#' gp <- create.gpData(pheno = dat, pedigree = ped)
#' A <- kin(gp, ret = "add")
#'
#' # assuming h2=sigma2u/(sigma2u+sigma2)=0.5
#' # no REML fit possible due to the limited number of observations
#' y <- c(132, 147, 156, 172)
#' names(y) <- paste("ID", 1:4, sep = "")
#' mod1 <- list(fit = list(sigma = c(1, 1), X = matrix(1, ncol = 1, nrow = 4)), kin = A, model = "BLUP", y = y, m = NULL)
#' # matrix A included all individuals (including those which should be predicted)
#' class(mod1) <- "gpMod"
#' predict(mod1, c("ID5", "ID6"))
#'
#' # prediction 'by hand'
#' X <- matrix(1, ncol = 1, nrow = 4)
#' Z <- diag(6)[-c(1, 2), ]
#' AI <- solve(A)
#' RI <- diag(4)
#'
#' res <- MME(X, Z, AI, RI, y)
#' res$u[1:2]
#' \dontrun{
#' # prediction of genetic performance of the last 50 individuals in the maize data set
#' data(maize)
#' maizeC <- codeGeno(maize)
#' U <- kin(maizeC, ret = "realized")
#' maizeC2 <- discard.individuals(maizeC, rownames(maizeC$pheno)[1201:1250])
#' modU <- gpMod(maizeC2, model = "BLUP", kin = U)
#' predict(modU, rownames(maizeC$pheno)[1201:1250])
#' }
#'
#' @export
#' @importFrom MASS ginv
#'
predict.gpMod <- function(object, newdata = NULL, ...) {
  if (class(object) != "gpMod") stop("'object' must be of class 'gpMod'")
  if (is.null(newdata)) {
    prediction <- object$g
  } else {
    model <- object$model

    if (model %in% c("BL", "BRR")) {
      if (!is.null(object$kin)) ("including a polygenic effect is not yet implemented")
      if (class(newdata) != "gpData") stop("object 'newdata' must be of class 'gpData'")
      X <- newdata$geno
      m <- object$markerEffects
      prediction <- X %*% m
    }
    if (model == "BLUP" & !is.null(object$markerEffects)) { # if marker effects are available
      if (class(newdata) != "gpData") stop("object 'newdata' must be of class 'gpData'")
      X <- newdata$geno
      m <- object$markerEffects
      prediction <- X %*% m
    }
    if (model == "BLUP" & is.null(object$markerEffects)) {
      #      prediction <- gpData$geno %*% t(gpData$geno[rownames(kin), ]) %*% ginv(kin) %*% genVal[rownames(kin)]
      #      prediction <- prediction[!names(prediction) %in% names(genVal)] / mean(prediction[names(genVal)]/genVal)
      if (any(newdata %in% names(object$y))) warning("Some individuals in newdata have been used also for model training")
      G <- object$kin[c(names(object$y), setdiff(y = names(object$y), x = newdata)), c(names(object$y), setdiff(y = names(object$y), x = newdata))]
      y <- object$y
      n <- length(y) # size of the training set
      Z <- cbind(diag(n), matrix(data = 0, nrow = n, ncol = length(setdiff(y = names(object$y), x = newdata))))
      X <- object$fit$X
      sigma2g <- object$fit$sigma[1]
      sigma2 <- object$fit$sigma[2]
      # diag(G) <- diag(G) + 0.00001
      GI <- ginv(G) * sigma2 / sigma2g # to ensure a solution for the inverse
      RI <- diag(n)
      sol <- MME(X, Z, GI, RI, y)
      prediction <- sol$u[-(1:n)]
      names(prediction) <- setdiff(y = names(object$y), x = newdata)
    }
  }
  return(prediction)
}

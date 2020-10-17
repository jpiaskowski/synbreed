# wrapper to apply genomic prediction models to an object of class gpData
# models used: regress/asreml and BLR

# date: 2011 - 01 - 16
# author: Valentin Wimmer
# date: 2011 - 05 - 03
# changes: Hans-Juergen Auinger
# date: 2011 - 11 - 21
# changes: return argument by Valentin Wimmer
# date: 2011 - 11 - 30
# changes: Hans-Juergen Auinger
# date: 2012 - 10 - 26
# changes: including R-asreml for BLUP by Hans-Juergen Auinger



#' Genomic predictions models for objects of class gpData
#'
#' This function fits genomic prediction models based on phenotypic and
#' genotypic data in an object of class \code{gpData}. The possible models are
#' Best Linear Unbiased Prediction (BLUP) using a pedigree-based or a
#' marker-based genetic relationship matrix and Bayesian Lasso (BL) or Bayesian
#' Ridge regression (BRR). BLUP models are fitted using the REML implementation
#' of the \code{regress} package (Clifford and McCullagh, 2012). The Bayesian
#' regression models are fitted using the Gibbs-Sampler of the \code{BGLR}
#' package (de los Campos and Perez, 2010). The covariance structure in the
#' BLUP model is defined by an object of class \code{relationshipMatrix}. The
#' training set for the model fit consists of all individuals with phenotypes
#' and genotypes. All data is restricted to individuals from the training set
#' used to fit the model.
#'
#' By default, an overall mean is added to the model. If no \code{kin} is
#' specified and \code{model = "BLUP"}, a G-BLUP model will be fitted. For
#' BLUP, further fixed and random effects can be added through the arguments
#' \code{fixed} and \code{random}.
#'
#' The marker effects \eqn{\hat{m}}{hatm} in the RR-BLUP model (available with
#' \code{markerEffects}) are calculated as \deqn{\hat{m}= X'G^{-1}\hat{g}}{m =
#' X'*Ginv*ghat} with \eqn{X} being the marker matrix, \eqn{G=XX'} and
#' \eqn{hat{g}}{ghat} the vector of predicted genetic values.
#'
#' Only a subset of the individuals - the training set - is used to fit the
#' model. This contains all individuals with phenotypes and genotypes. If
#' \code{kin} does not match the dimension of the training set (if, e.g.
#' ancestors are included), the respective rows and columns from the trainings
#' set are choosen.
#'
#' @param gpData object of class \code{gpData}
#' @param model \code{character}. Type of genomic prediction model.
#' \code{"BLUP"} indicates best linear unbiased prediction (BLUP) using REML
#' for both pedigree-based (P-BLUP) and marker-based (G-BLUP) model.
#' \code{"BL"} and \code{"BRR"} indicate Bayesian Lasso and Bayesian Ridge
#' Regression, respectively.
#' @param kin object of class \code{relationshipMatrix} (only required for
#' \code{model = "BLUP"}). Use a pedigree-based kinship to evaluate P-BLUP or a
#' marker-based kinship to evaluate G-BLUP. For \code{"BL"} and \code{"BRR"},
#' also a kinship structure may be used as additional polygenic effect \eqn{u}
#' in the Bayesian regression models (see \code{BGLR} package).
#' @param predict \code{logical}. If \code{TRUE}, genetic values will be
#' predicted for genotyped but not phenotyped individuals. Default is
#' \code{FALSE}. Note that this option is only meaningful for marker-based
#' models. For pedigree-based model, please use function \code{predict.gpMod}.
#' @param trait \code{numeric} or \code{character}. A vector with names or
#' numbers of the traits to fit the model
#' @param repl \code{numeric} or \code{character}. A vector with names or
#' numbers of the repeated values of \code{gpData$pheno} to fit the model
#' @param markerEffects \code{logical}.  Should marker effects be estimated for
#' a G-BLUP model, i.e. RR-BLUP? In this case, argument \code{kin} is ignored
#' (see Details). Plose note, that in this case also the variance components
#' pertaining to model G-BLUP are reported instead of those from the G-BLUP
#' model (see vignette). If the variance components are committed to
#' \code{crossVal}, it must be guaranteed that there also the RR-BLUP model is
#' used, e.g. no \code{cov.matrix} object should be specified.
#' @param fixed A formula for fixed effects. The details of model specification
#' are the same as for \code{lm} (only right hand side required). Only for
#' \code{model="BLUP"}.
#' @param random A formula for random effects of the model. Specifies the
#' matrices to include in the covariance structure. Each term is either a
#' symmetric matrix, or a factor. Independent Gaussian random effects are
#' included by passing the corresponding block factor. For mor details see
#' \code{\link[regress]{regress}}. Only for \code{model="BLUP"}
#' @param \dots further arguments to be used by the genomic prediction models,
#' i.e. prior values and MCMC options for the \code{BLR} function (see
#' \code{\link[BGLR]{BLR}}) or parameters for the REML algorithm in regress.
#'
#' @return Object of class \code{gpMod} which is a list of \item{fit}{The model
#' fit returned by the genomic prediction method} \item{model}{The model type,
#' see 'Arguments'} \item{y}{The phenotypic records for the individuals in the
#' training set} \item{g}{The predicted genetic values for the individuals in
#' the training set} \item{m}{Predicted SNP effects (if available)}
#' \item{kin}{Matrix \code{kin}}
#' @note The verbose output of the \code{BLR} function is written to a file
#' \code{BLRout.txt} in the working directory to prevent the screen output from
#' overload.
#' @author Valentin Wimmer, Hans-Juergen Auinger and Theresa Albrecht
#' @seealso \code{\link{kin}}, \code{\link{crossVal}}
#' @references Clifford D, McCullagh P (2012). regress: Gaussian Linear Models
#' with Linear Covariance Structure. R package version 1.3-8, URL
#' http://www.csiro.au.
#'
#' Gustavo de los Campos and Paulino Perez Rodriguez, (2010). BLR: Bayesian
#' Linear Regression. R package version 1.2.
#' http://CRAN.R-project.org/package=BGLR
#' @examples
#'
#' \dontrun{
#' library(synbreedData)
#' data(maize)
#' maizeC <- codeGeno(maize)
#'
#' # pedigree-based (expected) kinship matrix
#' K <- kin(maizeC, ret = "kin", DH = maize$covar$DH)
#'
#' # marker-based (realized) relationship matrix
#' # divide by an additional factor 2
#' # because for testcross prediction the kinship of DH lines is used
#' U <- kin(maizeC, ret = "realized") / 2
#' # BLUP models
#' # P-BLUP
#' mod1 <- gpMod(maizeC, model = "BLUP", kin = K)
#' # G-BLUP
#' mod2 <- gpMod(maizeC, model = "BLUP", kin = U)
#'
#' # Bayesian Lasso
#' prior <- list(varE = list(df = 3, S = 35), lambda = list(shape = 0.52, rate = 1e-4, value = 20, type = "random"))
#' mod3 <- gpMod(maizeC, model = "BL", prior = prior, nIter = 6000, burnIn = 1000, thin = 5)
#'
#' summary(mod1)
#' summary(mod2)
#' summary(mod3)
#' }
#'
#' @export gpMod
#' @import regress
#' @importFrom MASS ginv
#' @importFrom stats as.formula model.matrix na.pass terms
#' @importFrom utils capture.output
#'
gpMod <- function(gpData, model = c("BLUP", "BL", "BRR"), kin = NULL, predict = FALSE, trait = 1, repl = NULL, markerEffects = FALSE, fixed = NULL, random = NULL, ...) {

  # ASReml <- "package:asreml" %in% search()
  ans <- list()
  model <- match.arg(model)
  m <- NULL
  prediction <- NULL
  if (length(trait) > 1) cat("\n\tNOTE:   The return object will be a list of gpMod-objects!\n")
  if (is.null(fixed)) fixed <- ~1
  if (is.null(random)) {
    random <- " ~ "
    randomFormula <- ~1
  } else {
    randomFormula <- random
    random <- paste(paste(random, collapse = " "), "+ ")
  }
  if (model == "BLUP") {
    if (is.null(kin)) {
      if (!gpData$info$codeGeno) stop("Missing object 'kin', or use function codeGeno first!")
      # transposed crossproduct of the genotype matrix is used as relationship to obtain the variance components and mean of RR-BLUP
      if (markerEffects) kin <- tcrossprod(gpData$geno) else kin <- kin(gpData, ret = "realized")
    } else {
      if (markerEffects) kin <- tcrossprod(gpData$geno)
    }
  }
  for (i in trait) {
    df.trait <- gpData2data.frame(gpData, i, onlyPheno = TRUE, repl = repl)
    # take data from gpData object
    vec.bool <- colnames(df.trait) == "ID" | (colnames(df.trait) %in% dimnames(attr(terms(fixed), "factors"))[[1]]) | colnames(df.trait) %in% dimnames(attr(terms(randomFormula), "factors"))[[1]]
    cnt <- dim(gpData$pheno)[2]
    if (!is.null(gpData$phenoCovars)) cnt <- cnt + dim(gpData$phenoCovars)[2]
    if (i %in% 1:cnt) {
      yName <- dimnames(gpData$pheno)[[2]][as.numeric(i)]
      vec.bool[colnames(df.trait) %in% yName] <- TRUE
    } else {
      vec.bool <- vec.bool | colnames(df.trait) == i
      yName <- i
    }
    df.trait <- df.trait[, vec.bool]
    df.trait <- df.trait[!apply(is.na(df.trait), 1, sum), ]
    kinNames <- unique(df.trait$ID[!df.trait$ID %in% rownames(kin)])
    if (length(kinNames) != 0 & model == "BLUP") {
      df.trait <- df.trait[!df.trait$ID %in% kinNames, ]
      warning("Some phenotyped IDs are not in the kinship matrix!\nThese are removed from the analysis")
    }
    kinNew <- kin[unique(df.trait$ID), unique(df.trait$ID)]
    kinTS <- kin[df.trait$ID, df.trait$ID] # expand the matrix to what is needed
    if (model == "BLUP") {
      # if(!ASReml){
      res <- regress(as.formula(paste(yName, paste(fixed, collapse = " "))), Vformula = as.formula(paste(paste(random, collapse = " "), "kinTS")), data = df.trait, identity = TRUE, tol = 1e-8, ...)
      us <- BLUP(res)$Mean
      genVal <- us[grep("kinTS", names(us))]
      genVal <- genVal[!duplicated(names(genVal))]
      names(genVal) <- unlist(strsplit(names(genVal), "kinTS."))[(1:length(genVal)) * 2]
      # } else {
      #  kinTS <- kinNew[dimnames(gpData$pheno)[[1]], dimnames(gpData$pheno)[[1]]]
      #  covM.I <- try(solve(kinTS),TRUE)
      # adding constant to diagonal, if covM is singular
      #  if(class(covM.I)=="try-error"){
      #    warning("Covariance matrix is computationally singular: constant 1e-6 is added to the diagonal elements of the covariance matrix")
      #    covM.I <- solve(kinTS + diag(1e-6,ncol(kinTS)))
      # }
      # print warning in case of numerical problems
      # if(any(covM.I>1e8)) warning("Large >1e8 entries in the inverse covariance matrix")
      # kinGinv <- write.relationshipMatrix(covM.I,file=NULL,type="none",sorting="ASReml",digits=10)
      # attr(kinGinv, "rowNames") <- rownames(kinTS)

      # if("ginverse" %in% names(list(...)))
      #  asrObj <- asreml(as.formula(paste(yName, paste(fixed, collapse=""))), random=as.formula(paste(random, "giv(ID, var=TRUE)")), ginverse=c(list(ID=kinGinv), ginverse), data=df.trait)
      # else
      #  asrObj <- asreml(as.formula(paste(yName, paste(fixed, collapse=""))), random=as.formula(paste(random, "giv(ID, var=TRUE)")), ginverse=list(ID=kinGinv), data=df.trait)

      # genVal <-  asrObj$coefficients$random[substr(names(asrObj$coefficients$random), 1, 8) == "giv(ID)_"]
      # names(genVal) <- substr(names(genVal), 9, nchar(names(genVal)))
      # }
      # genVal <- NULL
      if (markerEffects) {
        m <- as.numeric(t(gpData$geno[rownames(kinNew), ]) %*% ginv(kinNew) %*% genVal[rownames(kinNew)])
        names(m) <- colnames(gpData$geno)
      }
      if (predict) {
        if (markerEffects) {
          prediction <- as.numeric(gpData$geno[!rownames(gpData$geno) %in% rownames(kinNew), ] %*% m)
          names(prediction) <- rownames(gpData$geno)[!rownames(gpData$geno) %in% rownames(kinNew)]
        } else {
          prediction <- gpData$geno %*% t(gpData$geno[rownames(kinNew), ]) %*% ginv(kinNew) %*% genVal[rownames(kinNew)]
          names(prediction) <- rownames(gpData$geno)
          prediction <- prediction[!dimnames(prediction)[[1]] %in% names(genVal)] / mean(prediction[names(genVal)] / genVal)
        }
      }
    }

    if (model == "BL") {
      if (random != " ~ ") stop("Random terms are not supported in model 'BL'!")
      if (dim(gpData$pheno)[3] > 1) stop("This method is not developed for a one-stage analysis yet. \nA phenotypic analysis have to be done fist.")
      X <- gpData$geno[rownames(gpData$geno) %in% df.trait$ID, ]
      y <- df.trait[df.trait$ID %in% rownames(gpData$geno), yName]
      if (fixed != " ~ 1") {
        XF <- as.formula(paste(yName, paste(fixed, collapse = " ")))
        XF <- model.matrix(XF, data = df.trait[df.trait$ID %in% rownames(gpData$geno), ], na.action = na.pass)
      } else {
        XF <- NULL
      }
      if (is.null(kin)) {
        if (is.null(XF)) {
          capture.output(res <- BLR(y = y, XL = X, ...), file = "BLRout.txt")
        } else {
          capture.output(res <- BLR(y = y, XL = X, XF = XF, ...), file = "BLRout.txt")
        }
      } else {
        if (is.null(XF)) {
          capture.output(res <- BLR(y = y, XL = X, GF = list(ID = df.trait$ID, A = kinTS), ...), file = "BLRout.txt")
        } else {
          capture.output(res <- BLR(y = y, XL = X, GF = list(ID = df.trait$ID, A = kinTS), XF = XF, ...), file = "BLRout.txt")
        }
      }
      genVal <- res$yHat - mean(res$yHat)
      names(genVal) <- rownames(X)
      m <- res$bL
      if (predict) {
        prediction <- as.numeric(gpData$geno[!rownames(gpData$geno) %in% names(genVal), ] %*% m)
        names(prediction) <- rownames(gpData$geno)[!rownames(gpData$geno) %in% names(genVal)]
      }
    }
    if (model == "BRR") {
      if (random != " ~ ") stop("Random terms are not supported in model 'BRR'!")
      if (dim(gpData$pheno)[3] > 1) stop("This method is not developed for a one-stage analysis yet. \nA phenotypic analysis has to be done fist.")
      X <- gpData$geno[rownames(gpData$geno) %in% df.trait$ID, ]
      y <- df.trait[df.trait$ID %in% rownames(gpData$geno), yName]
      if (fixed != " ~ 1") {
        XF <- as.formula(paste(yName, paste(fixed, collapse = " ")))
        XF <- model.matrix(XF, data = df.trait[df.trait$ID %in% rownames(gpData$geno), ], na.action = na.pass)
      } else {
        XF <- NULL
      }
      if (is.null(kin)) {
        if (is.null(XF)) {
          capture.output(res <- BLR(y = y, XR = X, ...), file = "BLRout.txt")
        } else {
          capture.output(res <- BLR(y = y, XR = X, XF = XF, ...), file = "BLRout.txt")
        }
      } else {
        if (is.null(XF)) {
          capture.output(res <- BLR(y = y, XR = X, GF = list(ID = df.trait$ID, A = kinTS), ...), file = "BLRout.txt")
        } else {
          capture.output(res <- BLR(y = y, XR = X, GF = list(ID = df.trait$ID, A = kinTS), XF = XF, ...), file = "BLRout.txt")
        }
      }
      genVal <- res$yHat - mean(res$yHat)
      names(genVal) <- rownames(X)
      m <- res$bR
      if (predict) {
        prediction <- as.numeric(gpData$geno[!rownames(gpData$geno) %in% names(genVal), ] %*% m)
        names(prediction) <- rownames(gpData$geno)[!rownames(gpData$geno) %in% names(genVal)]
      }
    }

    y <- df.trait[, yName]
    names(y) <- df.trait[, "ID"]
    ret <- list(fit = res, model = model, y = y, g = genVal, prediction = prediction, markerEffects = m, kin = kin)
    class(ret) <- "gpMod"
    ans[[i]] <- ret
    names(ans)[length(ans)] <- yName
  }
  if (length(trait) > 1) {
    class(ans) <- "gpModList"
    return(ans)
  } else {
    return(ret)
  }
}

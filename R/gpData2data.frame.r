# conversion of object class 'gpData' to data.frame



#' Merge of phenotypic and genotypic data
#'
#' Create a \code{data.frame} out of phenotypic and genotypic data in object of
#' class \code{gpData} by merging datasets using the common id. The shared data
#' set could either include individuals with phenotypes and genotypes (default)
#' or additional unphenotyped or ungenotyped individuals. In the latter cases,
#' the missing observations are filled by \code{NA}'s.
#'
#' Argument \code{all.geno} can be used to predict the genetic value of
#' individuals without phenotypic records using the \code{BGLR} package. Here,
#' the genetic value of individuals with \code{NA} as phenotype is predicted by
#' the marker profile.
#'
#' For multiple measures, phenotypic data in object \code{gpData} is arranged
#' with replicates in an \code{array}. With \code{gpData2data.frame} this could
#' be reshaped to "long" format with multiple observations in one column. In
#' this case, one column for the phenotype and 2 additional columns for the
#' \code{id} and the levels of the grouping variable (such as replications,
#' years of locations in multi-environment trials) are added.
#'
#' @param gpData object of class \code{gpData}
#' @param trait \code{numeric} or \code{character}. A vector with the names or
#' numbers of the trait that should be extracted from pheno. Default is
#' \code{1}.
#' @param onlyPheno scalar \code{logical}. Only return phenotypic data.
#' @param all.pheno scalar \code{logical}. Include all individuals with
#' phenotypes in the \code{data.frame} and fill the genotypic data with
#' \code{NA}.
#' @param all.geno scalar \code{logical}. Include all individuals with
#' genotypes in the \code{data.frame} and fill the phenotypic data with
#' \code{NA}.
#' @param repl \code{character} or \code{numeric}. A vector which contains
#' names or numbers of replication that should be drawn from the phenotypic
#' values and covariates. Default is \code{NULL}, i.e. all values are used.
#' @param phenoCovars \code{logical}. If \code{TRUE}, columns with the
#' phenotypic covariables are attached from element \code{phenoCovars} to the
#' \code{data.frame}. Only required for repeated measurements.
#' @param ...  further arguments to be used in function \code{reshape}. The
#' argument \code{times} could be useful to rename the levels of the grouping
#' variable (such as locations or environments).
#' @return A \code{data.frame} with the individuals names in the first column,
#' the phenotypes in the next column(s) and the marker genotypes in subsequent
#' columns.
#' @author Valentin Wimmer and Hans-Juergen Auinger
#' @seealso \code{\link{create.gpData}}, \code{\link{reshape}}
#' @keywords manip
#' @examples
#'
#' # example data with unrepeated observations
#' set.seed(311)
#'
#' # simulating genotypic and phenotypic data
#' pheno <- data.frame(Yield = rnorm(12, 100, 5), Height = rnorm(12, 100, 1))
#' rownames(pheno) <- letters[4:15]
#' geno <- matrix(sample(c("A", "A/B", "B", NA),
#'   size = 120, replace = TRUE,
#'   prob = c(0.6, 0.2, 0.1, 0.1)
#' ), nrow = 10)
#' rownames(geno) <- letters[1:10]
#' colnames(geno) <- paste("M", 1:12, sep = "")
#' # different subset of individuals in pheno and geno
#'
#' # create 'gpData' object
#' gp <- create.gpData(pheno = pheno, geno = geno)
#' summary(gp)
#' gp$covar
#'
#' # as data.frame with individuals with genotypes and phenotypes
#' gpData2data.frame(gp, trait = 1:2)
#' # as data.frame with all individuals with phenotypes
#' gpData2data.frame(gp, 1:2, all.pheno = TRUE)
#' # as data.frame with all individuals with genotypes
#' gpData2data.frame(gp, 1:2, all.geno = TRUE)
#'
#' # example with repeated observations
#' set.seed(311)
#'
#' # simulating genotypic and phenotypic data
#' pheno <- data.frame(ID = letters[1:10], Trait = c(
#'   rnorm(10, 1, 2), rnorm(10, 2, 0.2),
#'   rbeta(10, 2, 4)
#' ), repl = rep(1:3, each = 10))
#' geno <- matrix(rep(c(1, 0, 2), 10), nrow = 10)
#' colnames(geno) <- c("M1", "M2", "M3")
#' rownames(geno) <- letters[1:10]
#'
#' # create 'gpData' object
#' gp <- create.gpData(pheno = pheno, geno = geno, repeated = "repl")
#'
#' # reshape of phenotypic data and merge of genotypic data,
#' # levels of grouping variable loc are named "a", "b" and "c"
#' gpData2data.frame(gp, onlyPheno = FALSE, times = letters[1:3])
#' @export gpData2data.frame
#' @importFrom doBy orderBy
#' @importFrom methods as
#' @importFrom stats rnorm
#'
gpData2data.frame <- function(gpData, trait = 1, onlyPheno = FALSE, all.pheno = FALSE, all.geno = FALSE, repl = NULL, phenoCovars = TRUE, ...) {

  # check for class
  if (class(gpData) != "gpData") stop("object '", substitute(gpData), "' not of class 'gpData'")
  if (is.null(dimnames(gpData$pheno))) stop("There are not dimnames for the pheno slot of ", substitute(gpData))
  IDs <- dimnames(gpData$pheno)[[1]]
  reps <- dimnames(gpData$pheno)[[3]]
  pheno <- abind(gpData$pheno, matrix(1:dim(gpData$pheno)[1] + 10**ceiling(log10(dim(gpData$pheno)[1])), ncol = dim(gpData$pheno)[3], nrow = dim(gpData$pheno)[1], byrow = FALSE), along = 2)
  dimnames(pheno)[[2]][dim(pheno)[2]] <- "ID"
  # choose Traits
  if (all(is.numeric(trait))) trait <- (dimnames(pheno)[[2]])[trait]
  pheno <- array(pheno[, c("ID", trait), ], dim = c(dim(pheno)[1], length(trait) + 1, dim(pheno)[3]))
  dimnames(pheno) <- list(IDs, c("ID", trait), reps)
  # look for covariables
  if (phenoCovars) {
    if (!is.null(gpData$phenoCovars)) {
      pheno <- abind(pheno, gpData$phenoCovars, along = 2)
    }
  }
  # append column for each replication
  if (dim(pheno)[3] > 1) {
    pheno <- abind(matrix(rep(reps, each = dim(gpData$pheno)[1]), ncol = dim(gpData$pheno)[3], nrow = dim(gpData$pheno)[1], byrow = FALSE), pheno, along = 2)
    dimnames(pheno)[[2]][1] <- "repl"
  }
  if (!is.null(repl)) {
    if (all(repl %in% 1:dim(pheno)[3])) {
      repl <- dimnames(pheno)[[3]][repl]
    } else if (!all(repl %in% dimnames(pheno)[[3]])) stop("wrong replication names used")
  } else {
    repl <- dimnames(pheno)[[3]]
  }
  pheno <- as.data.frame(apply(pheno[, , repl], 2, cbind))
  if (phenoCovars) {
    for (i in names(gpData$info$attrPhenoCovars)) {
      if (gpData$info$attrPhenoCovars[i] == "numeric") {
        pheno[, i] <- as(as.character(pheno[, i]), gpData$info$attrPhenoCovars[i])
      }
    }
  }
  if (!is.null(repl)) {
    for (i in colnames(pheno)) {
      if (class(pheno[, i]) == "factor") {
        pheno[, i] <- as.factor(as.character(pheno[, i]))
      }
    }
  }
  pheno$ID <- IDs[as.numeric(as.factor(as.character(pheno$ID)))]
  if (!onlyPheno) {
    geno <- gpData$geno
    # merge genotypic and phenotypic data
    mergeData <- merge(pheno, geno, by.x = "ID", by.y = "row.names", all.x = all.pheno, all.y = all.geno)
    # omit row.names column
  } else {
    mergeData <- pheno
  }
  mergeData <- mergeData[, c(
    "ID", trait, ("repl")["repl" %in% colnames(mergeData)], colnames(mergeData)[colnames(mergeData) %in% colnames(gpData$geno)],
    colnames(mergeData)[colnames(mergeData) %in% dimnames(gpData$phenoCovars)[[2]]]
  )]
  # sort by ID
  for (i in trait) {
    mergeData[, i] <- as.numeric(as.character(mergeData[, i]))
  }
  if (all(mergeData$ID %in% 1:nrow(mergeData))) mergeData$ID <- as.numeric(mergeData$ID)
  if (is.null(mergeData$repl)) {
    mergeData <- orderBy(~ID, data = mergeData)
  } else {
    if (all(mergeData$repl %in% 1:dim(gpData$pheno)[3])) {
      mergeData$repl <- as.numeric(as.character(mergeData$repl))
    }
    mergeData <- orderBy(~ ID + repl, data = mergeData)
  }
  mergeData$ID <- as.character(mergeData$ID)
  rownames(mergeData) <- NULL
  return(mergeData)
}

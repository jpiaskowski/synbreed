# discard markers from class 'gpData'



#' Subsets for objects of class gpData
#'
#' The function produces subsets from an object of class \code{gpData} with
#' reduced markers. Marker informartion will be discarded from elements
#' \code{geno} and \code{map}
#'
#'
#' @param gpData object of class \code{gpData}
#' @param which character vector identifying names of markers which get
#' discarded in \code{geno} from a \code{gpData}-object.
#' @param whichNot character vector identifying names of markers which get kept
#' in \code{geno} from a \code{gpData}-object. Overwrites argument
#' \code{which}!
#' @return Object of class \code{gpData}
#' @author Valentin Wimmer and Hans-Juergen Auinger
#' @seealso \code{\link{create.gpData}}, \code{\link{add.markers}},
#' \code{\link{add.individuals}}, \code{\link{discard.individuals}}
#' @keywords manip
#' @examples
#'
#'
#' # example data
#' set.seed(311)
#' pheno <- data.frame(Yield = rnorm(10, 200, 5), Height = rnorm(10, 100, 1))
#' rownames(pheno) <- letters[1:10]
#' geno <- matrix(sample(c("A", "A/B", "B", NA),
#'   size = 120, replace = TRUE,
#'   prob = c(0.6, 0.2, 0.1, 0.1)
#' ), nrow = 10)
#' rownames(geno) <- letters[1:10]
#' colnames(geno) <- paste("M", 1:12, sep = "")
#' # one SNP is not mapped (M5)
#' map <- data.frame(chr = rep(1:3, each = 4), pos = rep(1:12))
#' map <- map[-5, ]
#' rownames(map) <- paste("M", c(1:4, 6:12), sep = "")
#' gp <- create.gpData(pheno = pheno, geno = geno, map = map)
#' summary(gp)
#'
#' # remove unmapped SNP M5 (which has no postion in the map)
#' gp2 <- discard.markers(gp, "M5")
#' summary(gp2)
#' \dontrun{
#' # add one new DH line to maize data
#' library(synbreedData)
#' data(maize)
#'
#' # delete markers
#' maize2 <- discard.individuals(maize, colnames(maize$geno)[1:50])
#' summary(maize2)
#' }
#'
#' @export discard.markers
#' @importFrom stats rnorm
#'
discard.markers <- function(gpData, which = NULL, whichNot = NULL) {
  if (class(gpData) != "gpData") {
    stop(substitute(gpData), " is not a gpData-object!")
  }
  # update geno
  if (is.null(whichNot)) {
    whichNot <- colnames(gpData$geno)[!colnames(gpData$geno) %in% which]
  }
  if (!is.null(attr(gpData$geno, "identical"))) attrG <- attr(gpData$geno, "identical") else attrG <- NULL
  gpData$geno <- gpData$geno[, colnames(gpData$geno) %in% whichNot]
  if (!is.null(attrG)) attr(gpData$geno, "identical") <- attrG
  # update map
  gpData$map <- gpData$map[rownames(gpData$map) %in% whichNot, ]
  return(gpData)
}

# discard individuals from class 'gpData'



#' Subsets for objects of class gpData
#'
#' The function produce subsets from an object of class \code{gpData} with
#' reduced individuals. Individual information will be discarded from elements
#' \code{geno}, \code{pheno}, \code{covar} and \code{pedigree}.
#'
#'
#' @param gpData object of class \code{gpData}
#' @param which character vector identifying names of individuals get discarded
#' from a \code{gpData}-object.
#' @param keepPedigree \code{logical}. Should the individual only be removed
#' from elements \code{geno} and \code{pheno} but kept in the pedigree?
#' @param whichNot character vector identifying names of individuals get kept
#' in a \code{gpData}-object. Overwrites argument \code{which}!
#' @return Object of class \code{gpData}
#' @author Valentin Wimmer and Hans-Juergen Auinger
#' @seealso \code{\link{create.gpData}}, \code{\link{add.individuals}},
#' \code{\link{add.markers}}, \code{\link{discard.markers}}
#' @keywords manip
#' @examples
#'
#'
#' # example data
#' set.seed(311)
#' pheno <- data.frame(Yield = rnorm(10, 200, 5), Height = rnorm(10, 100, 1))
#' rownames(pheno) <- letters[1:10]
#' geno <- matrix(sample(c("A", "A/B", "B", NA),
#'   size = 120, replace = TRUE,
#'   prob = c(0.6, 0.2, 0.1, 0.1)
#' ), nrow = 10)
#' rownames(geno) <- letters[1:10]
#' colnames(geno) <- paste("M", 1:12, sep = "")
#' # one SNP is not mapped (M5)
#' map <- data.frame(chr = rep(1:3, each = 4), pos = rep(1:12))
#' map <- map[-5, ]
#' rownames(map) <- paste("M", c(1:4, 6:12), sep = "")
#' gp <- create.gpData(pheno = pheno, geno = geno, map = map)
#' summary(gp)
#'
#' # discard genotypes with missing values in the marker matrix
#' gp3 <- discard.individuals(gp, names(which(rowSums(is.na(gp$geno)) > 0)))
#' summary(gp3)
#' \dontrun{
#' # add one new DH line to maize data
#' library(synbreedData)
#' data(maize)
#'
#' # delete individual
#' maize2 <- discard.individuals(maize, rownames(maize$geno)[1:10])
#' summary(maize2)
#' }
#'
#' @export discard.individuals
discard.individuals <- function(gpData, which = NULL, keepPedigree = FALSE, whichNot = NULL) {
  if (class(gpData) != "gpData") {
    stop(substitute(gpData), " is not a gpData-object!")
  }
  if (!is.null(whichNot)) {
    which <- c(
      rownames(gpData$geno)[!rownames(gpData$geno) %in% whichNot],
      dimnames(gpData$pheno)[[1]][[!dimnames(gpData$pheno)[[1]] %in% whichNot]]
    )
  }
  if (!all(which %in% gpData$covar$id)) stop("Some individuals are not in the ", substitute(gpData), "-object!")
  if (!is.null(gpData$pheno)) {
    # updata pheno
    phenoNames <- dimnames(gpData$pheno)
    phenoNames[[1]] <- phenoNames[[1]][!phenoNames[[1]] %in% which]
    phenoDim <- dim(gpData$pheno)
    phenoDim[1] <- length(phenoNames[[1]])
    gpData$pheno <- array(gpData$pheno[!rownames(gpData$pheno) %in% which, , ], dim = phenoDim)
    dimnames(gpData$pheno) <- phenoNames
  }
  # update geno
  if (!is.null(attr(gpData$geno, "identical"))) attrG <- attr(gpData$geno, "identical") else attrG <- NULL
  gpData$geno <- subset(gpData$geno, !rownames(gpData$geno) %in% which)
  if (!is.null(attrG)) attr(gpData$geno, "identical") <- attrG
  if (!is.null(gpData$phenoCovars)) {
    phenoCovarsNames <- dimnames(gpData$phenoCovars)
    phenoCovarsNames[[1]] <- phenoCovarsNames[[1]][!phenoCovarsNames[[1]] %in% which]
    phenoCovarsDim <- dim(gpData$phenoCovars)
    phenoCovarsDim[1] <- length(phenoCovarsNames[[1]])
    gpData$phenoCovars <- array(gpData$phenoCovars[!rownames(gpData$phenoCovars) %in% which, , ], dim = phenoCovarsDim)
    dimnames(gpData$phenoCovars) <- phenoCovarsNames
  }

  # update pedigree
  if (!is.null(gpData$pedigree)) {
    if (!keepPedigree) {
      gpData$pedigree <- subset(gpData$pedigree, !gpData$pedigree$ID %in% which)
      gpData$pedigree$Par1[gpData$pedigree$Par1 %in% which] <- 0
      gpData$pedigree$Par2[gpData$pedigree$Par2 %in% which] <- 0
      gpData$covar <- subset(gpData$covar, !gpData$covar$id %in% which)
    } else {
      gpData$covar[gpData$covar$id %in% which, c("phenotyped", "genotyped")] <- NA
    }
  }
  # update covar
  gpData$covar <- subset(gpData$covar, !gpData$covar$id %in% which)
  #  gpData$info$codeGeno <- FALSE

  return(gpData)
}

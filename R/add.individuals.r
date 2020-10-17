#' Add new individuals to objects of class gpData
#'
#' This function extends an object of class \code{gpData} by adding new
#' phenotypes, genotypes and pedigree.
#'
#' \code{colnames} in \code{geno}, \code{pheno} and \code{pedigree} must match
#' existing names in \code{gpData} object.
#'
#' @param gpData object of class \code{gpData} to be updated
#' @param pheno \code{data.frame} with new rows for phenotypes with
#' \code{rownames} indicating individuals. For repeated values the ID should be
#' stored in a column with name \code{"ID"}.
#' @param geno \code{matrix} with new rows for genotypic data with
#' \code{rownames} indicating individuals
#' @param pedigree \code{data.frame} with new rows for pedigree data
#' @param covar \code{data.frame} with new rows for \code{covar} information
#' with \code{rownames} indicating individuals
#' @param repl The column of the pheno \code{data.frame} for the replicated
#' measures. If the values are not repeated or this column is named
#' \code{"repl"} this argument is not needed.
#' @return object of class \code{gpData} with new individuals
#' @author Valentin Wimmer
#' @seealso \code{\link{add.markers}}, \code{\link{discard.individuals}}
#' @keywords manip
#' @examples
#'
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
#' # new phenotypic data
#' newPheno <- data.frame(Yield = 200, Height = 100, row.names = "newLine")
#' # simulating genotypic data
#' newGeno <- matrix(sample(c("A", "A/B", "B"), ncol(gp$geno), replace = TRUE), nrow = 1)
#' rownames(newGeno) <- "newLine"
#' # new pedigree
#' newPedigree <- create.pedigree(ID = "newLine", Par1 = 0, Par2 = 0, gener = 0)
#'
#' gp2 <- add.individuals(gp, pheno = newPheno, geno = newGeno, pedigree = newPedigree)
#' \dontrun{
#' # add one new DH line to maize data
#' library(synbreedData)
#' data(maize)
#' newDHpheno <- data.frame(Trait = 1000, row.names = "newDH")
#' # simulating genotypic data
#' newDHgeno <- matrix(sample(c(0, 1), ncol(maize$geno), replace = TRUE), nrow = 1)
#' rownames(newDHgeno) <- "newDH"
#' # new pedigree
#' newDHpedigree <- create.pedigree(ID = "newDH", Par1 = 0, Par2 = 0, gener = 0)
#' # new covar information
#' newDHcovar <- data.frame(family = NA, DH = 1, tbv = 1000, row.names = "newDH")
#'
#' # add individual
#' maize2 <- add.individuals(maize, newDHpheno, newDHgeno, newDHpedigree, newDHcovar)
#' summary(maize2)
#' }
#'
#' @export add.individuals
#' @importFrom doBy orderBy
#' @importFrom stats rnorm
#'
add.individuals <- function(gpData, pheno = NULL, geno = NULL, pedigree = NULL, covar = NULL, repl = NULL) {
  # check if individuals are allready in data
  if (any(rownames(pheno) %in% gpData$covar$id) | any(rownames(geno) %in% gpData$covar$id) |
    any(pheno$ID %in% gpData$covar$id)) {
    stop("some of the individuals of are already in ", substitute(gpData))
  }
  if (is.null(repl)) {
    repl <- "repl"
  } else if (!all(unique(pheno[, repl]) %in% dimnames(gpData$pheno)[[3]])) stop("Your values for replication is not in the dimnames of ", substitute(gpData$pheno))
  colnames(pheno)[colnames(pheno) == repl] <- "repl"
  # merge phenotypic data
  if (dim(gpData$pheno)[3] == 1) {
    if (!"ID" %in% colnames(pheno)) pheno$ID <- rownames(pheno)
    repl <- NULL
  } else {
    if (!"ID" %in% colnames(pheno)) stop("In", substitute(pheno), "the columns 'ID' and/or 'repl' are/is missing!")
  }
  if (!is.null(pheno)) {
    if (!all(!colnames(pheno)[colnames(pheno) != "ID"] %in% dimnames(gpData$pheno)[[2]] |
      !colnames(pheno)[colnames(pheno) != "ID"] %in% dimnames(gpData$phenoCovars)[[2]])) {
      stop("different phenotypes (colnames) in '", substitute(gpData$pheno), "' or '", substitute(gpData$phenoCovar), "' and '", substitute(pheno), "'")
    }
  }
  df.pheno <- gpData2data.frame(gpData, onlyPheno = TRUE, trait = dimnames(gpData$pheno)[[2]])
  if (!all(colnames(df.pheno) %in% colnames(pheno))) warning("Not all traits and covariates are available in the new data!")
  pheno[, colnames(df.pheno)[!colnames(df.pheno) %in% colnames(pheno)]] <- NA
  pheno <- pheno[, colnames(df.pheno)]
  pheno <- rbind(df.pheno, pheno)
  rm(df.pheno)
  if (dim(gpData$pheno)[3] == 1) {
    rownames(pheno) <- pheno$ID
    pheno$ID <- NULL
  }
  # merge genotypic data
  if (!is.null(geno)) if (any(colnames(geno) != colnames(gpData$geno))) stop("different markers (colnames) in 'gpData$geno' and 'geno'")
  geno <- rbind(gpData$geno, geno)

  # merge pedigree
  pedigree <- rbind(gpData$pedigree, pedigree)
  # reorder if necessary
  if (!is.null(pedigree)) pedigree <- orderBy(~gener, pedigree)

  # merge covar  (not with first columns)

  if (any(colnames(covar) %in% c("genotyped", "phenotyped", "id"))) stop("not specify columns 'genotyped','phenotyped' and 'id' in 'covar' ")
  if (!is.null(covar)) {
    cc <- data.frame(gpData$covar[, !colnames(gpData$covar) %in% c("genotyped", "phenotyped", "id")])
    rownames(cc) <- rownames(gpData$covar)
    colnames(cc) <- colnames(gpData$covar)[!colnames(gpData$covar) %in% c("genotyped", "phenotyped", "id")]
    covarUpdate <- rbind(cc, covar)
  }
  else {
    covarUpdate <- gpData$covar
  }

  # need id in rownames for create.gpData
  rownames(covarUpdate) <- c(as.character(gpData$covar$id), rownames(covar))

  # create new gpData object
  ret <- create.gpData(pheno = pheno, geno = geno, map = gpData$map, pedigree = pedigree, covar = covarUpdate, map.unit = gpData$info$map.unit, modCovar = dimnames(gpData$phenoCovars)[[2]], repeated = repl)
  return(ret)
}

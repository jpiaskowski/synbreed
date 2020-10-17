#' Conversion between objects of class 'cross' and 'gpData'
#'
#' Function to convert an object of class \code{gpData} to an object of class
#' \code{cross} (F2 intercross class in the package \code{qtl}) and vice versa.
#' If not done before, function \code{codeGeno} is used for recoding in
#' \code{gpData2cross}.
#'
#' In \code{cross}, genotypic data is splitted into chromosomes while in
#' \code{gpData} genotypic data comprises all chromosomes because separation
#' into chromosomes in not required for genomic prediction. Note that coding of
#' genotypic data differs between classes. In \code{gpData}, genotypic data is
#' coded as the number of copies of the minor allele, i.e. 0, 1 and 2. Thus,
#' function \code{codeGeno} should be applied to \code{gpData} before using
#' \code{gpData2cross} to ensure correct coding. In \code{cross}, coding for F2
#' intercross is: AA = 1, AB = 2, BB = 3. When using \code{gpData2cross} or
#' \code{cross2gpData}, resulting genotypic data has correct format.
#'
#' @aliases gpData2cross cross2gpData
#' @param gpData object of class \code{gpData} with non-empty elements for
#' \code{pheno}, \code{geno} and \code{map}
#' @param ...  further arguments for function \code{codeGeno}. Only used in
#' \code{gpData2cross}.
#' @return Object of class \code{cross} of \code{gpData} for function
#' \code{gpData2cross} and \code{cross2gpData}, respectively.
#' @author Valentin Wimmer and Hans-Juergen Auinger
#' @seealso \code{\link{create.gpData}}, \code{\link[qtl]{read.cross}} ,
#' \code{\link{codeGeno}}
#' @references Broman, K. W. and Churchill, S. S. (2003). R/qtl: Qtl mapping in
#' experimental crosses. Bioinformatics, (19):889-890.
#' @examples
#'
#' \dontrun{
#' library(synbreedData)
#' # from gpData to cross
#' data(maize)
#' maizeC <- codeGeno(maize)
#' maize.cross <- gpData2cross(maizeC)
#' # descriptive statistics
#' summary(maize.cross)
#' plot(maize.cross)
#'
#' # use function scanone
#' maize.cross <- calc.genoprob(maize.cross, step = 2.5)
#' result <- scanone(maize.cross, pheno.col = 1, method = "em")
#' # display of LOD curve along the chromosome
#' plot(result)
#'
#'
#' # from cross to gpData
#' data(fake.f2)
#' fake.f2.gpData <- cross2gpData(fake.f2)
#' summary(fake.f2.gpData)
#' }
#'
#' @export gpData2cross
#' @importFrom qtl read.cross
#'
gpData2cross <- function(gpData, ...) {
  # check for class
  if (class(gpData) != "gpData") stop("object '", substitute(gpData), "' not of class 'gpData'")
  # check on geno and map
  if (is.null(gpData$geno) | is.null(gpData$map)) stop("'geno' and 'map' needed in", substitute(gpData))
  if (dim(gpData$pheno)[3] > 1) {
    stop("You can only use unreplicated values for cross!")
  } else {
    # use codeGeno if not yet done
    if (!gpData$info$codeGeno) stop("Use function codeGeno before gpData2cross!")

    # only use individuals with genotypes and phenotypes
    genoPheno <- gpData$covar$id[gpData$covar$genotyped & gpData$covar$phenotyped]
    # read information from gpData
    geno <- data.frame(gpData$geno[rownames(gpData$geno) %in% genoPheno, ])
    phenoDim <- dim(gpData$pheno)
    phenoNames <- dimnames(gpData$pheno)
    phenoDim[1] <- sum(dimnames(gpData$pheno)[[1]] %in% genoPheno)
    phenoNames[[1]] <- dimnames(gpData$pheno)[[1]][dimnames(gpData$pheno)[[1]] %in% genoPheno]
    pheno <- gpData$pheno[dimnames(gpData$pheno)[[1]] %in% genoPheno, , ]
    pheno <- array(pheno, dim = phenoDim)
    dimnames(pheno) <- phenoNames
    pheno <- apply(pheno, 2, rbind) # possible because of unreplicated data!!!
    rownames(pheno) <- phenoNames[[1]]
    if (dim(gpData$pheno)[3] > 1) pheno$repl <- rep(1:dim(gpData$pheno)[3], each = dim(gpData$pheno)[1])
    if (!is.null(gpData$phenoCovars)) pheno <- cbind(pheno, data.frame(apply(gpData$phenoCovars[dimnames(gpData$phenoCovars)[[1]] %in% genoPheno, , ], 2, rbind)))
    map <- gpData$map
    n <- nrow(geno)
    pheno <- as.data.frame(pheno)
  }

  # split markers (+pos) and genotypes on chromosomes
  genoList <- split(cbind(rownames(map), map$pos, t(geno)), map$chr)
  # result is a list
  # function to bring each list element in right format
  addData <- function(x) {
    ret <- list()
    Nm <- length(x) / (n + 2)

    # elements of x:
    #  1:Nm: marker names
    #  (Nm+1):(2*Nm): marker positions
    #  rest: genotypes as vector

    # add 1 to genotypes
    # coding for F2 intercross: AA=1, AB=2, BB=3
    ret[["data"]] <- matrix(as.numeric(x[-(1:(2 * Nm))]) + 1, nrow = n, ncol = Nm, byrow = TRUE, dimnames = list(NULL, x[1:Nm]))
    ret[["map"]] <- as.numeric(x[(Nm + 1):(2 * Nm)])
    names(ret[["map"]]) <- x[1:Nm]
    # this may have to be changed
    class(ret) <- "A"
    ret
  }
  #
  # apply function to each list element
  genoList <- lapply(genoList, addData)
  # create object 'cross'
  cross <- list(geno = genoList, pheno = pheno)
  class(cross) <- c("f2", "cross")
  cross
}

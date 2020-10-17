#' Prepare genotypic data for Beagle
#'
#' Create input file for Beagle software (Browning and Browning 2009) from an
#' object of class \code{gpData}. This function is created for usage within
#' function \code{codeGeno} to impute missing values.
#'
#' The Beagle software must be used chromosomewise. Consequently, \code{gp}
#' should contain only data from one chromosome (use \code{discard.markers, see
#' Examples}).
#'
#' @param gp \code{gpData} object with elements \code{geno} and \code{map}
#' @param wdir \code{character}. Directory for Beagle input files
#' @param prefix \code{character}. Prefix for Beagle input files.
#' @return No value is returned. Function creates files
#' \code{[prefix]ingput.bgl} with genotypic data in Beagle input format and
#' \code{[prefix]marker.txt} with marker information used by Beagle.
#' @author Valentin Wimmer
#' @seealso \code{\link{codeGeno}}
#' @references B L Browning and S R Browning (2009) A unified approach to
#' genotype imputation and haplotype phase inference for large data sets of
#' trios and unrelated individuals. Am J Hum Genet 84:210-22
#' @keywords manip
#' @examples
#'
#' map <- data.frame(chr = c(1, 1, 1, 1, 1, 2, 2, 2, 2), pos = 1:9)
#' geno <- matrix(sample(c(0, 1, 2, NA), size = 10 * 9, replace = TRUE), nrow = 10, ncol = 9)
#' colnames(geno) <- rownames(map) <- paste("SNP", 1:9, sep = "")
#' rownames(geno) <- paste("ID", 1:10 + 100, sep = "")
#'
#' gp <- create.gpData(geno = geno, map = map)
#' gp1 <- discard.markers(gp, rownames(map[map$chr != 1, ]))
#' \dontrun{
#' write.beagle(gp1, prefix = "test")
#' }
#'
#' @export write.beagle
#' @importFrom utils str write.table
#'
write.beagle <- function(gp, wdir = getwd(), prefix) {

  # information from arguments
  geno <- gp$geno
  n <- nrow(geno)
  M <- ncol(geno)
  if (is.null(n) | is.null(M)) {
    print(str(geno))
    stop("Wrong genotypic information!")
  }

  # prepare input file input.bgl
  if (any(grep(" ", rownames(geno)))) stop("no blanks allowed in rownames(geno) when running beagle")
  firstLine <- c("I", "id", rep(rownames(geno), each = 2))
  # assume coding AA, AB, BB
  allele1 <- substr(geno, 1, 1)
  allele2 <- substr(geno, 2, 2)
  hap <- rep(NA, 2 * n)
  hap[seq(1, (2 * n * M), by = 2)] <- allele1
  hap[seq(2, (2 * n * M), by = 2)] <- allele2
  # markers = rows
  haplotypes <- matrix(hap, nrow = M, ncol = 2 * n, byrow = TRUE)
  # write input file
  input.bgl <- rbind(firstLine, cbind("M", colnames(geno), haplotypes))
  write.table(input.bgl, file = file.path(wdir, paste(prefix, "input.bgl", sep = "")), quote = FALSE, col.names = FALSE, row.names = FALSE)

  # write marker file
  getAlleles <- function(Genotype) paste(unique(unlist(strsplit(paste(as.character(Genotype[!is.na(Genotype)])), split = ""))), collapse = " ")
  write.table(cbind(rownames(gp$map), gp$map$pos, apply(geno, 2, getAlleles)), file = file.path(wdir, paste(prefix, "marker.txt", sep = "")), quote = FALSE, col.names = FALSE, row.names = FALSE)
}

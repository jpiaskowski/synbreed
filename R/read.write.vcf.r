#' Prepare genotypic data in vcf-Format
#'
#' Create vcf-file for miscellaneous applications. Within the package it is
#' used to write files for beagle usage.
#'
#' The function writes a vcf-file. The format of the output is "GT". Other
#' formats are not supported.
#'
#' @param gp \code{gpData} object with elements \code{geno} and \code{map}
#' @param file \code{character}. Filename for writing the file.
#' @param unphased \code{logical}. The default is TRUE. Than the seperator
#' between the alleles is \code{"/"}, and the possible codings are \code{"0/0"}
#' for \code{0} in the genotype matrix, \code{"0/1"} for \code{1} and
#' \code{"1/1"} for 2. For getting a phased output, use \code{unphased=FALSE}.
#' Than the seperator is \code{"|"}. For hetercygous genotypes you have to
#' change the 1 to -1 if you like to get the coding \code{"1|0"}, So posible
#' codings in this case are \code{"0|0"} for \code{0} in the genotype matrix,
#' \code{"0|1"} for \code{1}, \code{"1|0"} for \code{-1} and \code{"1|1"} for
#' 2.
#' @return No value is returned. Function creates files
#' \code{[prefix]ingput.bgl} with genotypic data in Beagle input format and
#' \code{[prefix]marker.txt} with marker information used by Beagle.
#' @author Hans-Juergen Auinger
#' @seealso \code{\link{read.vcf2matrix}}, \code{\link{read.vcf2list}},
#' \code{\link{codeGeno}}
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
#' write.vcf(gp1, file = "test.vcf")
#' }
#'
#' @export write.vcf
write.vcf <- function(gp, file, unphased = TRUE) {
  if (is.null(nrow(gp$geno)) | is.null(ncol(gp$geno))) {
    print(str(gp$geno))
    stop("Wrong genotypic information!")
  }
  if (unphased) sepSign <- "/" else sepSign <- "|"
  if (gp$info$map.unit == "kb") gp$map$pos <- gp$map$pos * 1000
  if (gp$info$map.unit == "Mb") gp$map$pos <- gp$map$pos * 1000000
  if (!gp$info$map.unit %in% c("bp", "kb", "Mb")) {
    gp$map$pos <- "."
    warning("Positions will not be written to the file! Only basepair positions will be printed")
  } else if (any((gp$map$pos - round(gp$map$pos, digits = 0)) > 1e-6)) stop("Your map positions and the map.unit do not fit!")
  geno <- as.data.frame(t(gp$geno), stringsAsFactors = FALSE)
  s00 <- paste("0", "0", sep = sepSign)
  s01 <- paste("0", "1", sep = sepSign)
  s11 <- paste("1", "1", sep = sepSign)
  s10 <- paste("1", "0", sep = sepSign)
  geno[geno == 0] <- s00
  geno[geno == 1] <- s01
  geno[geno == 2] <- s11
  geno[geno == -1] <- s10
  bgl <- cbind(
    data.frame(CHROM = paste("chr", gp$map$chr, sep = ""), POS = formatC(gp$map$pos, format = "f", digits = 0), ID = rownames(gp$map), REF = "A", ALT = "G", QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT", stringsAsFactors = FALSE),
    geno
  )
  if (any(grep(" ", colnames(geno)))) stop("no blanks allowed in IDs!")
  cat(file = file, "##fileformat=VCFv4.1\n##filedate=", Sys.Date(), '\n##source="write.vcf of R-synbreed"\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#')
  cat(file = file, paste(colnames(bgl), collapse = "\t"), "\n", append = TRUE)
  write.table(bgl,
    file = file,
    quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, sep = "\t", na = paste(".", ".", sep = sepSign)
  )
}



#' Read data of a vcf-file to a matrix
#'
#' To easily read genomic data in vcf-Format to a matrix. Function
#' \code{codeGeno} uses \code{read.vcf2matrix} with imputing by beagle.
#'
#'
#' @param file \code{character}. The name of the file which the data are to be
#' read from.
#' @param FORMAT \code{character}. The default is \code{"GT"}. If there are
#' more formats in your vcf-file you can decide which one you like to have in
#' your output matrix.
#' @param coding This option has only an effect with \code{FORMAT="GT"}.
#' \code{allele} gives you back the alles as defined as REF and ALT in your
#' vcf-file. \code{ref} gives you back \code{"0"} for the reference allele and
#' \code{"1"} for the alternative allele.
#' @param IDinRow \code{logical}. Default is \code{TRUE}, this means the
#' genotypes are in the rows and the markers in the column. For \code{FALSE} it
#' is the other way round.
#' @param cores \code{numeric}. Specifies the number of cores for parallel
#' computing.
#' @return A matrix (\code{\link[base]{matrix}}) containing a representation of
#' the data in the file.
#' @author Hans-Juergen Auinger
#' @seealso \code{\link{write.vcf}}, \code{\link{read.vcf2list}}
#' @examples
#'
#' \dontrun{
#' library(synbreedData)
#' data(maize)
#' maize$info$map.unit <- "kb"
#' maize <- codeGeno(maize)
#' write.vcf(maize, "maize.vcf")
#' geno <- read.vcf2matrix("maize.vcf")
#' }
#'
#' @export read.vcf2matrix
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster parLapply stopCluster mclapply
#' @importFrom utils read.table str write.table
#'
read.vcf2matrix <- function(file, FORMAT = "GT", coding = c("allele", "ref"), IDinRow = TRUE, cores = 1) {
  multiLapply <- function(x, y, ..., cores = cores) {
    if (.Platform$OS.type == "windows" & cores > 1) {
      cl <- makeCluster(min(cores, detectCores()))
      registerDoParallel(cl)
      parLapply(cl, x, y, ...)
      stopCluster(cl)
    } else {
      mclapply(x, y, ..., mc.cores = cores)
    }
  }
  coding <- match.arg(coding)
  cnt <- 0
  while (scan(file = file, what = "character", skip = cnt, nlines = 1, quiet = TRUE)[1] != "#CHROM") cnt <- cnt + 1
  Mnames <- scan(file, what = "character", skip = cnt, nlines = 1, quiet = TRUE)
  geno <- read.table(file = file, sep = "\t", header = FALSE, skip = cnt + 1, stringsAsFactors = FALSE)
  colnames(geno) <- Mnames
  rownames(geno) <- geno$ID
  ref <- geno$REF
  alternative <- geno$ALT
  form <- unlist(strsplit(geno$FORMAT, ":"))
  geno <- geno[, !colnames(geno) %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")]
  geno[1:nrow(geno), 1:ncol(geno)] <- unlist(multiLapply(geno, strsplit, ":", cores = cores))[rep(form, ncol(geno)) == FORMAT]
  if (FORMAT == "GT" & coding == "allele") {
    geno[geno == "0|0"] <- rep(paste(ref, ref, sep = "|"), ncol(geno))[geno == "0|0"]
    geno[geno == "1|0"] <- rep(paste(alternative, ref, sep = "|"), ncol(geno))[geno == "1|0"]
    geno[geno == "0|1"] <- rep(paste(ref, alternative, sep = "|"), ncol(geno))[geno == "0|1"]
    geno[geno == "1|1"] <- rep(paste(alternative, alternative, sep = "|"), ncol(geno))[geno == "1|1"]
  }
  if (IDinRow) geno <- t(geno)
  return(geno)
}



#' Read data of a vcf-file to a matrix
#'
#' Function for easily read genomic data in vcf-Format to a list, which
#' contains the map information and the marker information.
#'
#'
#' @param file \code{character}. The name of the file which the data are to be
#' read from.
#' @param FORMAT \code{character}. The default is \code{"GT"}. If there are
#' more formats in your vcf-file you can decide which one you like to have in
#' your output matrix.
#' @param coding This option has only an effect with \code{FORMAT="GT"}.
#' \code{allele} gives you back the alles as defined as REF and ALT in your
#' vcf-file. \code{ref} gives you back \code{"0"} for the reference allele and
#' \code{"1"} for the alternative allele.
#' @param IDinRow \code{logical}. Default is \code{TRUE}, this means the
#' genotypes are in the rows and the markers in the column. For \code{FALSE} it
#' is the other way round.
#' @param cores \code{numeric}. Specifies the number of cores for parallel
#' computing.
#' @return A list with a matrix (\code{\link[base]{matrix}}) containing a
#' representation of the genotypic data in the file and a map of classes GenMap
#' and data.frame.
#' @author Hans-Juergen Auinger
#' @seealso \code{\link{write.vcf}}, \code{\link{read.vcf2matrix}}
#' @examples
#'
#' \dontrun{
#' library(synbreedData)
#' data(maize)
#' maize$info$map.unit <- "kb"
#' maize <- codeGeno(maize)
#' write.vcf(maize, "maize.vcf")
#' genInfo <- read.vcf2list("maize.vcf")
#' }
#'
#' @export read.vcf2list
read.vcf2list <- function(file, FORMAT = "GT", coding = c("allele", "ref"), IDinRow = TRUE, cores = 1) {
  multiLapply <- function(x, y, ..., cores = cores) {
    if (.Platform$OS.type == "windows" & cores > 1) {
      cl <- makeCluster(min(cores, detectCores()))
      registerDoParallel(cl)
      parLapply(cl, x, y, ...)
      stopCluster(cl)
    } else {
      mclapply(x, y, ..., mc.cores = cores)
    }
  }
  coding <- match.arg(coding)
  cnt <- 0
  while (scan(file = file, what = "character", skip = cnt, nlines = 1, quiet = TRUE)[1] != "#CHROM") cnt <- cnt + 1
  Mnames <- scan(file, what = "character", skip = cnt, nlines = 1, quiet = TRUE)
  geno <- read.table(file = file, sep = "\t", header = FALSE, skip = cnt + 1, stringsAsFactors = FALSE)
  colnames(geno) <- Mnames
  rownames(geno) <- geno$ID
  ref <- geno$REF
  alternative <- geno$ALT
  form <- unlist(strsplit(geno$FORMAT, ":"))
  map <- geno[, c("#CHROM", "POS")]
  colnames(map) <- c("chr", "pos")
  class(map) <- c("GenMap", "data.frame")
  geno <- geno[, !colnames(geno) %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")]
  geno[1:nrow(geno), 1:ncol(geno)] <- unlist(multiLapply(geno, strsplit, ":", cores = cores))[rep(form, ncol(geno)) == FORMAT]
  if (FORMAT == "GT" & coding == "allele") {
    geno[geno == "0|0"] <- rep(paste(ref, ref, sep = "|"), ncol(geno))[geno == "0|0"]
    geno[geno == "1|0"] <- rep(paste(alternative, ref, sep = "|"), ncol(geno))[geno == "1|0"]
    geno[geno == "0|1"] <- rep(paste(ref, alternative, sep = "|"), ncol(geno))[geno == "0|1"]
    geno[geno == "1|1"] <- rep(paste(alternative, alternative, sep = "|"), ncol(geno))[geno == "1|1"]
  }
  if (IDinRow) geno <- t(geno)
  return(list(geno = geno, map = map))
}

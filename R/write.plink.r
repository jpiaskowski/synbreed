#' Prepare data for PLINK
#'
#' Create input files and corresponding script for PLINK (Purcell et al. 2007)
#' to estimate pairwise LD through function \code{pairwiseLD}.
#'
#'
#' @param gp \code{gpData} object with elements \code{geno} and \code{map}
#' @param wdir \code{character}. Directory for PLINK input files
#' @param prefix \code{character}. Prefix for PLINK input files.
#' @param ld.threshold \code{numeric}. Threshold for the LD used in PLINK.
#' @param type \code{character}. Specifies the type of return value for PLINK.
#' @param ld.window \code{numeric}. Window size for pairwise differences which
#' will be reported by PLINK (only for \code{use.plink=TRUE}; argument
#' \code{--ld-window-kb} in PLINK) to thin the output dimensions. Only SNP
#' pairs with a distance < \code{ld.window} are reported (default = 99999).
#' @return No value returned. Files \code{prefix.map}, \code{prefix.ped} and
#' \code{prefixPlinkScript.txt} are created in the working directory
#' @author Valentin Wimmer
#' @seealso \code{\link{pairwiseLD}}
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender
#' D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) PLINK: a
#' toolset for whole-genome association and population-based linkage analysis.
#' American Journal of Human Genetics, 81.
#' @keywords manip
#' @examples
#'
#' \dontrun{
#' library(synbreedData)
#' write.plink(maize, type = "data.frame")
#' }
#'
#' @export write.plink
#' @importFrom utils write.table
#'
write.plink <- function(gp, wdir = getwd(), prefix = paste(substitute(gp)), ld.threshold = 0, type = c("data.frame", "matrix"), ld.window = 99999) {

  # information from arguments
  geno <- gp$geno
  n <- nrow(geno)
  M <- ncol(geno)

  # recode genotypic data
  geno[geno == 0] <- "BB" ## ? or "AA"
  geno[geno == 2] <- "AA"
  geno[geno == 1] <- "AB"

  # prepare .ped file
  # family ID, individual ID, paternal ID, maternal ID, sex, phenotype, genotypedata
  # we omit cols 3-6
  if (is.null(family)) {
    warning("missing family information in gpData$covar")
    gp$covar$family <- 1
  }

  family <- gp$covar$family[gp$covar$genotyped] # only for individuals with genotypes

  id <- gp$covar$id[gp$covar$genotyped] # only for individuals with genotypes
  # combine in one file
  ped <- data.frame(family, id, geno)

  # prepare .map file
  # Chromosome ID, SNP ID, genetic distance, bp pos
  chr <- gp$map$chr
  if (length(chr) != M) stop("PLINK does not allow different number of markers in 'map' and 'geno'")
  snpID <- rownames(gp$map)
  gPos <- bpPos <- gp$map$pos
  # else gPos <- rep(0,M)
  # if(gp$info$map.unit=="bp") bpPos <- gp$map$pos
  # else bpPos <- rep(0,M)
  # combine in one file
  map <- data.frame(chr, snpID, gPos, bpPos)

  # write data for PLINK
  write.table(ped, file = file.path(wdir, paste(prefix, ".ped", sep = "")), quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(map, file = file.path(wdir, paste(prefix, ".map", sep = "")), quote = FALSE, col.names = FALSE, row.names = FALSE)

  # write PLINK script (sink does not work for Linux)
  # zz <- file(file.path(wdir,paste(prefix,"plinkScript.txt",sep="")), open="wt")
  # sink(zz)
  # cat("--ped ",prefix,".ped \n",sep="",append=TRUE)
  # cat("--map ",prefix,".map \n",sep="",append=TRUE)
  # cat("--compound-genotypes \n")
  # cat("--out ",prefix,"\n")
  # cat("--no-parents \n")
  # cat("--no-sex \n")
  # cat("--no-pheno \n")
  # cat("--allow-no-sex \n")
  # if(type=="matrix") cat("--matrix \n")
  # else {
  # cat("--ld-window 99999\n")
  # cat("--ld-window-r2 ",ld.threshold)
  # }
  # cat("--r2 \n")

  # sink()
  # unlink(file.path(wdir,paste(prefix,"plinkScript.txt",sep="")))

  cat("--ped ", prefix, ".ped \n", "--map ", prefix, ".map \n", "--compound-genotypes \n", "--out ", prefix, "\n", "--no-parents \n", "--no-sex \n", "--no-pheno \n", "--allow-no-sex \n", ifelse(type == "matrix", "--matrix \n", paste("--ld-window ", nrow(map), "\n", "--ld-window-kb ", ld.window, "\n", "--ld-window-r2 ", ld.threshold, "\n", sep = "")), "--r2 \n", sep = "", file = file.path(wdir, paste(prefix, "plinkScript.txt", sep = "")))
}

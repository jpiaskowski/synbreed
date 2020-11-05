
# library(synbreedData)
# data(maize)
# gpData1 <- maize
# gpData2 <- maize
# gpData2$covar$id <- as.character(gpData2$covar$id)
# gpData2$covar$id[gpData2$covar$id %in% rownames(gpData2$geno)] <- paste(gpData2$covar$id[gpData2$covar$id %in% rownames(gpData2$geno)], 0, sep="")
# gpData2$pedigree$ID[gpData2$pedigree$ID %in% rownames(gpData2$geno)] <- paste(gpData2$pedigree$ID[gpData2$pedigree$ID %in% rownames(gpData2$geno)], 0, sep="")
# rownames(gpData2$geno) <- paste(rownames(gpData2$geno), 0, sep="")
# dimnames(gpData2$pheno)[[1]] <- paste(dimnames(gpData2$pheno)[[1]], 0, sep="")
# gpData2$pheno <- abind(gpData2$pheno, gpData2$pheno, along =2)
# dimnames(gpData2$pheno)[[2]] <- c("IDTrait", "Trait2")


#' Join two \code{gpData} objects
#'
#' Function for joining two \code{gpData} objects
#'
#' The function writes a vcf file. The format of the output is "GT". Other
#' formats are not supported.
#'
#' @param gpData1 A \code{gpData} object with at least elements \code{geno} and
#' \code{map}
#' @param gpData2 Second \code{gpData} object with at least elements
#' \code{geno} and \code{map}
#' @return A \code{gpData} object is returned containing \code{gpData1}
#' and \code{gpData2}
#' @author Hans-Juergen Auinger
#' @seealso \code{\link{create.gpData}} \code{\link{codeGeno}}
#' @examples
#'
#' \dontrun{
#' add.gpData(maize, maize)
#' }
#'
#' @export add.gpData
#' @importFrom abind abind
#'
add.gpData <- function(gpData1, gpData2) {
  if (!is.null(gpData1$info$version)) stop(paste("Recode ", substitute(gpData1), "! You have used an old version to create/code ", substitute(gpData1), sep = ""))
  if (substr(gpData1$info$version, 47, 50) < 0.12) stop(paste("Recode ", substitute(gpData1), "! You have used an old version to create/code ", substitute(gpData1), sep = ""))
  if (!is.null(gpData2$info$version)) stop(paste("Recode ", substitute(gpData1), "! You have used an old version to create/code ", substitute(gpData1), sep = ""))
  if (substr(gpData2$info$version, 47, 50) < 0.12) stop(paste("Recode ", substitute(gpData2), "! You have used an old version to create/code ", substitute(gpData2), sep = ""))
  if (is.null(gpData1$map)) {
    if (is.null(gpData2$map)) {
      map <- NULL
    } else {
      map <- gpData2$map
    }
  } else {
    if (is.null(gpData2$map)) {
      geno <- gpData1$map
    } else {
      if (gpData1$info$map.unit != gpData2$info$map.unit) stop("map.units in gpData-objects are different!")
      map1 <- gpData1$map[!is.na(gpData1$map$chr), ]
      map2 <- gpData2$map[!is.na(gpData2$map$chr), ]
      map1$names <- rownames(map1)
      map2$names <- rownames(map2)
      map <- unique(rbind(map1, map2))
      if (any(duplicated(map$names))) stop("Marker map information are different between the gpData-objects!")
      map$names <- NULL
    }
  }
  if (!is.null(gpData1$pheno)) {
    pheno1 <- gpData2data.frame(gpData1, onlyPheno = TRUE, trait = 1:dim(gpData1$pheno)[2], stringsAsFactors = TRUE)
  } else {
    pheno1 <- NULL
  }
  if (!is.null(gpData1$pheno)) {
    pheno2 <- gpData2data.frame(gpData2, onlyPheno = TRUE, trait = 1:dim(gpData2$pheno)[2], stringsAsFactors = TRUE)
  } else {
    pheno2 <- NULL
  }
  if (is.null(pheno1)) {
    if (is.null(pheno2)) {
      pheno <- NULL
    } else {
      pheno <- pheno2
    }
  } else {
    if (is.null(pheno2)) {
      pheno <- pheno1
    } else {
      if (ncol(pheno1) == ncol(pheno2)) cnmsSwtch <- !all(colnames(pheno1) == colnames(pheno2))
      if (ncol(pheno1) != ncol(pheno2) | cnmsSwtch) {
        colPhen1 <- colnames(pheno2)[!colnames(pheno2) %in% colnames(pheno1)]
        colPhen2 <- colnames(pheno1)[!colnames(pheno1) %in% colnames(pheno2)]
        if (!is.null(colPhen1)) {
          pheno1[, colPhen1] <- NA
        }
        if (!is.null(colPhen2)) {
          pheno2[, colPhen2] <- NA
        }
        pheno2 <- pheno2[, colnames(pheno1)]
      }
      pheno <- rbind(pheno1, pheno2)
    }
  }
  if (is.null(gpData1$geno)) {
    if (is.null(gpData2$geno)) {
      geno <- NULL
    } else {
      geno <- gpData2$geno
    }
  } else {
    if (is.null(gpData2$geno)) {
      geno <- gpData1$geno
    } else {
      if (ncol(gpData1$geno) == ncol(gpData2$geno)) cnmsSwtch <- !all(colnames(gpData1$geno) == colnames(gpData2$geno))
      if (ncol(gpData1$geno) != ncol(gpData2$geno) | cnmsSwtch) {
        colGen1 <- colnames(gpData2$geno)[!colnames(gpData2$geno) %in% colnames(gpData1$geno)]
        colGen2 <- colnames(gpData1$geno)[!colnames(gpData1$geno) %in% colnames(gpData2$geno)]
        if (!is.null(colGen1)) {
          gpData1$geno[, colGen1] <- NA
        }
        if (!is.null(colGen2)) {
          gpData2$geno[, colGen2] <- NA
        }
        gpData2$geno <- gpData2$geno[, colnames(gpData1$geno)]
      }
      geno <- rbind(gpData1$geno, gpData2$geno)
    }
  }
  if (is.null(gpData1$pedigree)) {
    if (is.null(gpData2$pedigree)) {
      pedigree <- NULL
    } else {
      pedigree <- gpData2$pedigree
    }
  } else {
    if (is.null(gpData2$pedigree)) {
      pedigree <- gpData1$pedigree
    } else {
      pedigree <- add.pedigree(gpData1$pedigree, gpData2$pedigree)
    }
  }
  create.gpData(
    geno = geno, pheno = pheno, map = map, pedigree = pedigree, # covar=covarUpdate,
    map.unit = gpData1$info$map.unit, modCovar = dimnames(gpData1$phenoCovars)[[2]] # ,repeated=repl
  )
  return(0)
}

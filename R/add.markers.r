#' Add new markers to an object of class gpData
#'
#' This function adds new markers to the element \code{geno} of an object of
#' class \code{gpData} and updates the marker map.
#'
#' \code{rownames} in argument \code{geno} must match \code{rownames} in the
#' element \code{geno} object of class \code{gpData}.
#'
#' @param gpData object of class \code{gpData} to be updated
#' @param geno \code{matrix} with new columns
#' @param map \code{data.frame} with columns 'chr' and 'pos' for new markers
#' @return object of class \code{gpData} with new markers
#' @author Valentin Wimmer
#' @seealso \code{\link{add.individuals}}, \code{\link{discard.markers}}
#' @keywords manip
#' @examples
#'
#' # creating gpData object
#' # phenotypic data
#' pheno <- data.frame(Yield = rnorm(10, 100, 5), Height = rnorm(10, 10, 1))
#' rownames(pheno) <- 1:10
#' # genotypic data
#' geno <- matrix(sample(c(1, 0, 2, NA),
#'   size = 120, replace = TRUE,
#'   prob = c(0.6, 0.2, 0.1, 0.1)
#' ), nrow = 10)
#' rownames(geno) <- 1:10
#' # genetic map
#' map <- data.frame(chr = rep(1:3, each = 4), pos = rep(1:12))
#' colnames(geno) <- rownames(map) <- paste("M", 1:12, sep = "")
#' # as gpData object
#' gp <- create.gpData(pheno, geno, map)
#'
#'
#' # new data
#' geno2 <- matrix(c(0, 0, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 0, 2, 1, 1, 1, 2, 2, 2), ncol = 2)
#' rownames(geno2) <- 1:10
#' map2 <- data.frame(pos = c(0.3, 5), chr = c(1, 2))
#' rownames(map2) <- colnames(geno2) <- c("M13", "M14")
#'
#' # adding new markers
#' gp2 <- add.markers(gp, geno2, map2)
#' summary(gp2)
#' summary(gp)
#' @export add.markers
#' @importFrom doBy orderBy
#' @importFrom stats rnorm
#'
add.markers <- function(gpData, geno, map = NULL) {
  class(gpData$map) <- "data.frame"
  # check if markers are allready in data
  if (any(colnames(geno) %in% colnames(gpData$geno)) | any(rownames(map) %in% colnames(gpData$geno))) {
    stop("some of the markers of ", substitute(geno), " are allready in ", substitute(gpData))
  }
  if (is.null(gpData$map) & !is.null(map)) stop("There is no map available for ", substitute(gpData), "!")
  if (!all(rownames(geno) %in% gpData$covar$id)) stop("You like to put new individuals into the data set!\nUse add.individuals() for that!")
  if (!all(rownames(map) %in% colnames(geno))) stop("There are markers in the map, which don't have information in ", substitute(geno), "!")
  # take names form map if available
  if (is.null(colnames(geno)) & !is.null(map)) {
    if (ncol(geno) == nrow(map)) {
      colnames(geno) <- rownames(map)
    } else {
      stop("Check the colnames of", substitute(geno), " and the rownames of ", substitute(map), "!")
    }
  }

  # merge genotypic data
  nmiss <- sum(!rownames(gpData$geno) %in% rownames(geno))
  if (nmiss > 0) {
    geno <- rbind(matrix(NA, nrow = nmiss, ncol = ncol(geno)), geno)
    rownames(geno)[1:nmiss] <- rownames(gpData$geno)[!rownames(gpData$geno) %in% rownames(geno)]
  }
  geno <- geno[match(rownames(gpData$geno), rownames(geno)), ]
  gpData$geno <- cbind(gpData$geno, geno)
  # first column as rownames and delete first column
  # merge map
  if (is.null(map)) {
    # map <- gpData$map[1:ncol(geno),]
    map <- data.frame(chr = rep(NA, ncol(geno)), pos = rep(NA, ncol(geno)))
    rownames(map) <- colnames(geno)
  } else if (nrow(map) != ncol(geno)) {
    map[colnames(geno)[!colnames(geno) %in% rownames(map)], ] <- NA
  }
  map <- rbind(gpData$map, map)
  # same approach as in create.gpData
  map$sor <- substr(map$chr, nchar(as.character(map$chr)), nchar(as.character(map$chr)))
  if (!all(!unique(map$sor)[!is.na(unique(map$sor))] %in% 0:9)) map$sor <- 1
  map <- map[order(as.character(rownames(map))), ]
  map <- orderBy(~ sor + chr + pos, data = map)
  map$sor <- NULL
  gpData$map <- map
  gpData$geno <- gpData$geno[, match(rownames(gpData$map), colnames(gpData$geno))]
  gpData$info$codeGeno <- FALSE
  # sortcolumns in geno, too
  gpData$geno <- gpData$geno[, rownames(gpData$map)]
  # create new gpData object
  class(gpData$map) <- "GenMap"
  return(gpData)
}

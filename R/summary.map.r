#' Summary of marker map information
#'
#' This function can be used to summarize information from a marker map in an
#' object of class \code{gpData}.  Return value is a \code{data.frame} with one
#' row for each chromosome and one row summarizing all chromosomes.
#'
#' Summary statistics of differences are based on euclidian distances between
#' markers with non-missing position in \code{map}, i.e. \code{pos!=NA}.
#'
#' @param map \code{data.frame} with columns \code{chr} and \code{pos} or a
#' \code{gpData} object with element \code{map}
#' @param cores \code{numeric}. Specifies the number of cores for parallel
#' computing.
#' @return A \code{data.frame} with one row for each chromosome and the
#' intersection of all chromosomes and columns \item{noM}{number of markers}
#' \item{range}{range of positions, i.e. difference between first and last
#' marker} \item{avDist}{avarage distance of markers} \item{maxDist}{maximum
#' distance of markers} \item{minDist}{minimum distance of markers}
#' @author Valentin Wimmer
#' @seealso \code{\link{create.gpData}}
#' @examples
#'
#' \dontrun{
#' library(synbreedData)
#' data(maize)
#' summaryGenMap(maize)
#' }
#'
#' @export summaryGenMap
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster parLapply stopCluster mclapply
#' @importFrom stats weighted.mean
#'
summaryGenMap <- function(map, cores = 1) {
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
  # information from arguments
  if (class(map) == "gpData") map <- map$map
  if (is.null(map)) stop("No map available")
  chr <- map$chr
  pos <- map$pos
  # extract information
  # number of markers
  len <- tapply(pos, chr, length)
  rge <- tapply(pos, chr, max, na.rm = TRUE) - tapply(pos, chr, min, na.rm = TRUE)
  # differences of markers on each chr
  diffs <- tapply(pos[!is.na(pos)], chr[!is.na(pos)], diff, na.rm = TRUE)
  avDist <- as.numeric(multiLapply(diffs, mean, na.rm = TRUE, cores = cores))
  maxDist <- as.numeric(multiLapply(diffs, max, na.rm = TRUE, cores = cores))
  minDist <- as.numeric(multiLapply(diffs, min, na.rm = TRUE, cores = cores))
  # return data.frame
  ret <- data.frame(noM = len, length = rge, avDist = avDist, maxDist = maxDist, minDist = minDist, row.names = names(len))
  # keep same order as original map
  ret <- ret[order(order(unique(chr))), ]
  # sum over all chr
  all <- data.frame(noM = sum(ret$noM, na.rm = TRUE), length = sum(ret$length, na.rm = TRUE), avDist = weighted.mean(ret$avDist, ret$noM, na.rm = TRUE), maxDist = max(ret$maxDist, na.rm = TRUE), minDist = min(ret$minDist, na.rm = TRUE))
  rownames(all) <- paste(rownames(ret)[1], "-", rownames(ret)[nrow(ret)])
  ret <- rbind(ret, all)
  return(ret)
}

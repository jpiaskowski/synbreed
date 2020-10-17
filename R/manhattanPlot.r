#' Manhattan plot for SNP effects
#'
#' Plot of SNP effects along the chromosome, e.g. for the visualization of
#' marker effects generated by function \code{gpMod}.
#'
#'
#' @param b object of class \code{gpMod} with marker effects or numeric vector
#' of marker effects to plot
#' @param gpData object of class \code{gpData} with map position
#' @param colored \code{logical}. Color the chromosomes?. The \code{default} is
#' \code{FALSE} with chromosomes distinguished by grey tones.
#' @param add If \code{TRUE}, the plot is added to an existing plot. The
#' default is \code{FALSE}.
#' @param pch a vector of plotting characters or symbols: see
#' \code{\link[graphics]{points}}. The default is an open circle.
#' @param ylab a title for the y axis: see \code{\link[graphics]{title}}.
#' @param colP color for the different chromosomes and points.
#' @param \dots further arguments for function \code{plot}
#' @author Valentin Wimmer and Hans-Juergen Auinger
#' @keywords hplot
#' @examples
#'
#' \dontrun{
#' library(synbreedData)
#' data(mice)
#' # plot only random noise
#' b <- rexp(ncol(mice$geno), 3)
#' manhattanPlot(b, mice)
#' }
#'
#' @export manhattanPlot
#' @importFrom grDevices grey rainbow
#' @importFrom methods is
#' @importFrom graphics axis box plot points
#'
manhattanPlot <- function(b, gpData = NULL, colored = FALSE, add = FALSE, pch = 19, ylab = NULL, colP = NULL, ...) {
  if (is.null(gpData)) {
    plot(b, ...)
  } else {
    if (is.null(gpData$map)) stop("missing map in gpData object ", substitute(gpData))
    if (class(b) == "gpMod") b <- b$markerEffects
    b <- b[!(is.na(gpData$map$pos) | is.na(gpData$map$chr))]
    gpData$map <- gpData$map[!(is.na(gpData$map$pos) | is.na(gpData$map$chr)), ]
    if (is.null(colP)) {
      if (colored) {
        colP <- rainbow(6)
      } else {
        colP <- rep(c(grey(0.3), grey(0.7)), times = length(unique(gpData$map$chr)))
      }
      colP <- colP[(as.numeric(gpData$map$chr) - 1) %% 6 + 1]
    }
    chrs <- cumsum(tapply(gpData$map$pos, gpData$map$chr, max))
    namChrs <- names(chrs)
    chrs <- c(0, chrs[1:(length(chrs) - 1)])
    names(chrs) <- namChrs
    chr <- as.numeric(chrs[gpData$map$chr]) + as.numeric(gpData$map$pos) + as.numeric(as.factor(gpData$map$chr)) * 0.01
    if (!add) {
      plot(chr, b, col = colP, type = "p", axes = FALSE, pch = pch, ylab = ylab, ...)
      axis(side = 1, at = c(chr[!duplicated(gpData$map$chr)], max(chr, na.rm = TRUE)), labels = NA, cex = .9, lwd.ticks = 1, ...)
      axis(
        side = 1, at = chr[!duplicated(gpData$map$chr)] + diff(c(chr[!duplicated(gpData$map$chr)], max(chr, na.rm = TRUE))) / 2,
        tick = FALSE, labels = unique(gpData$map$chr), hadj = 0, padj = 0, ...
      )
      axis(side = 2, ...)
      box()
    } else {
      points(chr, b, col = colP[(as.numeric(gpData$map$chr) - 1) %% 6 + 1], pch = pch, ...)
    }
  }
}

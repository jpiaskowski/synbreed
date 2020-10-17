#' Visualization of pedigree
#'
#' A function to visualize pedigree structure by a graph using the
#' \code{igraph} package. Each genotype is represented as vertex and direct
#' offsprings are linked by an edge.
#'
#' The pedigree is structured top to bottom. The first generation is printed in
#' the first line. Links over more than one generation are possible as well as
#' genotypes with only one (known) parent. Usually, no structure in one
#' generation is plotted. If an \code{effect} is given, the genotypes are
#' ordered by this effect in the horizontal direction and a labeled axis is
#' plotted at the bottom.
#'
#' @aliases plot.pedigree
#'
#' @param x object of class \code{pedigree} or object of class \code{gpData}
#' with element \code{pedigree}
#' @param effect vector of length \code{nrow(pedigree)} with effects to plot on
#' the x axis
#' @param ...  Other arguments for function \code{igraph.plotting}
#' @return A named graph visualizing the pedigree structure. Color is used to
#' distinguish sex.
#' @note This function uses the plotting method for graphs in the library
#' \code{igraph}
#' @author Valentin Wimmer and Hans-Juergen Auinger
#' @seealso \code{\link{create.pedigree}}, \code{\link{simul.pedigree}}
#' @keywords hplot
#' @examples
#'
#' id <- paste("ID", 1:9, sep = "0")
#' par1 <- paste("ID", c("", "", "", "", 1, 1, 1, 4, 7), sep = "0")
#' par2 <- paste("ID", c("", "", "", "", 2, 3, 2, 5, 8), sep = "0")
#' ped1 <- create.pedigree(id, par1, par2, unknown = "ID0")
#' ped1
#' plot(ped1)
#'
#' # create 2nd pedigree object
#' Id <- paste("ID", 10:16, sep = "")
#' Par1 <- paste("ID", c("", "", 1, 1, 6, 7, 7), sep = "0")
#' Par2 <- paste("ID", c("", "", 10, "08", "09", 11, 14), sep = "")
#' ped2 <- create.pedigree(Id, Par1, Par2, unknown = c("ID0", "ID"))
#' ped2
#'
#' ped <- add.pedigree(ped1, ped2)
#' plot(ped)
#' @export
#' @importFrom igraph graph.data.frame
#' @importFrom methods is
#' @importFrom graphics axis plot
#'
plot.pedigree <- function(x, effect = NULL, ...) {

  # get data
  if (any(class(x) == "gpData")) pedigree <- x$pedigree else pedigree <- x

  # catch errors
  if (min(pedigree$gener) != 0) pedigree$gener <- pedigree$gener - min(pedigree$gener)
  if (!is.null(effect) & length(effect) != nrow(pedigree)) stop("length of effect does not equal nrow(pedigree)")



  # set parents which are not in pedigree to unknown
  pedigree[!(pedigree$Par1 %in% pedigree$ID), ]$Par1 <- 0
  pedigree[!(pedigree$Par2 %in% pedigree$ID), ]$Par2 <- 0

  # remove duplicated entries (same Par1 and Par2) for onlyFamily plot
  # if (onlyFamily) pedigree <- pedigree[!(duplicated(pedigree[,c("Par1","Par2")]) & pedigree$gener==max(pedigree$gener)) ,]

  # build pedigree graph
  relations <- rbind(as.matrix(pedigree[pedigree$Par1 != 0, c("Par1", "ID")], ncol = 2), as.matrix(pedigree[pedigree$Par2 != 0, c("Par2", "ID")], ncol = 2))
  ped.graph <- graph.data.frame(relations, directed = TRUE, pedigree)

  # extract generation
  gener <- pedigree$gener

  # define x-y positions
  n <- nrow(pedigree)
  pos <- matrix(data = NA, nrow = n, ncol = 2)
  pos[, 2] <- max(gener) - gener
  if (is.null(effect)) {
    pos[, 1] <- order(gener, partial = order(pedigree$ID, decreasing = TRUE)) - cumsum(c(0, table(gener)))[gener + 1]
  } else {
    pos[, 1] <- effect
  }

  # nice format
  myscale <- function(x) {
    if (length(x) == 1) {
      x <- 0
    } else {
      x <- unlist(scale(x))
    }
    return(x)
  }

  if (is.null(effect)) pos[n:1, 1] <- unlist(tapply(pos[, 1], pos[, 2], myscale))

  # color for vertex
  cols <- rep("lightblue", n)
  if (!is.null(pedigree$sex)) cols[pedigree$sex == 0] <- "palevioletred1"
  # plot
  # ps.options(family="serif")
  # if(onlyFamily) plot(ped.graph,rescale=TRUE,vertex.label=NA,layout=pos,edge.color=1,edge.width=0.5,edge.arrow.size=0,vertex.size=0,...)

  # else{
  plot(ped.graph, rescale = TRUE, vertex.label = pedigree$ID, layout = pos, edge.color = 1, edge.width = 0.5, edge.arrow.size = 0.5, vertex.label.family = "sans", vertex.color = cols, ...)
  # }
  if (!is.null(effect)) axis(side = 1, at = seq(-1, 1, length = 10), labels = round(seq(min(pos[, 1]), max(pos[, 1]), length = 10), 0))
  return(ped.graph)
}

#' Writing relationshipMatrix in table format
#'
#' This function can be used to write an object of class "relationshipMatrix"
#' in the table format used by other software, i.e. WOMBAT or ASReml.  The
#' resulting table has three columns, the row, the column and the entry of the
#' (inverse) relationshipMatrix.
#'
#' Note that "WOMBAT" only uses the generalized inverse relationship matrix and
#' expects a file with the name "ranef.gin", where 'ranef' is the name of the
#' random effect with option 'GIN' in the 'MODEL' part of the parameter file.
#' For ASREML, either the relationship could be saved as "*.grm" or its
#' generalized inverse as "*.giv".
#'
#' @param x Object of class "relationshipMatrix"
#' @param file Path where the output should be written . If \code{NULL} the
#' result is returned in R.
#' @param sorting Type of sorting. Use "WOMBAT" for 'row-wise' sorting of the
#' table and "ASReml" for 'column-wise' sorting.
#' @param type A character string indicating which form of
#' \code{relationshipMatrix} should be returned. One of "ginv" (Moore-Penrose
#' generalized inverse), "inv" (inverse), or "none" (no inverse).
#' @param digits Numeric.  The result is rounded to \code{digits}.
#' @author Valentin Wimmer
#' @references Meyer, K. (2006) WOMBAT - A tool for mixed model analyses in
#' quantitative genetics by REML, J. Zhejinag Uni SCIENCE B 8: 815-821.
#'
#' Gilmour, A., Cullis B., Welham S., and Thompson R. (2000) ASREML. program
#' user manual. NSW Agriculture, Orange Agricultural Institute, Forest Road,
#' Orange, Australia .
#' @keywords IO
#' @examples
#'
#' \dontrun{
#' # example with 9 individuals
#' id <- 1:9
#' par1 <- c(0, 0, 0, 0, 1, 1, 1, 4, 7)
#' par2 <- c(0, 0, 0, 0, 2, 3, 2, 5, 8)
#' gener <- c(0, 0, 0, 0, 1, 1, 1, 2, 3)
#' ped <- create.pedigree(id, par1, par2, gener)
#' gp <- create.gpData(pedigree = ped)
#'
#' A <- kin(ped, ret = "add")
#' write.relationshipMatrix(A, type = "ginv")
#' }
#'
#' @export write.relationshipMatrix
#' @importFrom MASS ginv
#' @importFrom utils write.table
#'
write.relationshipMatrix <- function(x, file = NULL, sorting = c("WOMBAT", "ASReml"), type = c("ginv", "inv", "none"), digits = 10) {
  type <- match.arg(type)
  sorting <- match.arg(sorting)

  if (sorting == "WOMBAT" & type != "ginv") stop("'type' must be 'ginv' for WOMBAT")

  # pass (inverse) relationship matrix
  if (type == "ginv") rMinv <- ginv(x)
  if (type == "inv") rMinv <- solve(x)
  if (type == "none") rMinv <- x

  rMinv <- round(rMinv, digits)

  # add rownames and colnames
  res <- data.frame(
    Row = rep(1:nrow(rMinv), nrow(rMinv)),
    Column = rep(1:nrow(rMinv), each = nrow(rMinv)),
    coeff = as.numeric(rMinv),
    lower = as.logical(lower.tri(rMinv, diag = TRUE))
  )
  rm(rMinv)



  # only use lower triangle
  res <- res[res$lower == TRUE, c("Row", "Column", "coeff")]

  if (sorting == "ASReml") {
    res <- res[order(res$Row, res$Column), ]
  }
  if (sorting == "WOMBAT") {
    res <- res[, c(2, 1, 3)]
    res <- res[order(res$Column, res$Row), ]
  }
  res <- res[res$coeff != 0, ]

  # write to given file
  if (!is.null(file)) {
    write.table(res, file, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
    rm(x)
  } else {
    attr(res, "rowNames") <- rownames(x)
    rm(x)
    return(res)
  }
}

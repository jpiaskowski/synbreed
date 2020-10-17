

#' Extract or replace part of map data.frame
#' 
#' Extract or replace part of an object of class \code{GenMap}.
#' 
#' 
#' @name [.GenMap
#' @docType data
#' @param x object of class ''GenMap''
#' @param ...  indices
#' @examples
#' 
#' \dontrun{
#' data(maize)
#' head(maize$map)
#' }
#' 
NULL





#' Extract or replace part of relationship matrix
#' 
#' Extract or replace part of an object of class \code{relationshipMatrix}.
#' 
#' 
#' @name [.relationshipMatrix
#' @docType data
#' @param x object of class ''relationshipMatrix''
#' @param ...  indices
#' @examples
#' 
#' \dontrun{
#' data(maize)
#' U <- kin(codeGeno(maize),ret="realized")
#' U[1:3,1:3]
#' }
#' 
NULL




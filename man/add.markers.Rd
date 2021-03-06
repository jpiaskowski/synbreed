% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.markers.r
\name{add.markers}
\alias{add.markers}
\title{Add new markers to an object of class gpData}
\usage{
add.markers(gpData, geno, map = NULL)
}
\arguments{
\item{gpData}{object of class \code{gpData} to be updated}

\item{geno}{\code{matrix} with new columns}

\item{map}{\code{data.frame} with columns 'chr' and 'pos' for new markers}
}
\value{
object of class \code{gpData} with new markers
}
\description{
This function adds new markers to the element \code{geno} of an object of
class \code{gpData} and updates the marker map.
}
\details{
\code{rownames} in argument \code{geno} must match \code{rownames} in the
element \code{geno} object of class \code{gpData}.
}
\examples{

# creating gpData object
# phenotypic data
pheno <- data.frame(Yield = rnorm(10, 100, 5), Height = rnorm(10, 10, 1))
rownames(pheno) <- 1:10
# genotypic data
geno <- matrix(sample(c(1, 0, 2, NA),
  size = 120, replace = TRUE,
  prob = c(0.6, 0.2, 0.1, 0.1)
), nrow = 10)
rownames(geno) <- 1:10
# genetic map
map <- data.frame(chr = rep(1:3, each = 4), pos = rep(1:12))
colnames(geno) <- rownames(map) <- paste("M", 1:12, sep = "")
# as gpData object
gp <- create.gpData(pheno, geno, map)


# new data
geno2 <- matrix(c(0, 0, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 0, 2, 1, 1, 1, 2, 2, 2), ncol = 2)
rownames(geno2) <- 1:10
map2 <- data.frame(pos = c(0.3, 5), chr = c(1, 2))
rownames(map2) <- colnames(geno2) <- c("M13", "M14")

# adding new markers
gp2 <- add.markers(gp, geno2, map2)
summary(gp2)
summary(gp)
}
\seealso{
\code{\link{add.individuals}}, \code{\link{discard.markers}}
}
\author{
Valentin Wimmer
}
\keyword{manip}

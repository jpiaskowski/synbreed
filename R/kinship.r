#' Relatedness based on pedigree or marker data
#'
#' This function implements different measures of relatedness between
#' individuals in an object of class \code{gpData}: (1) Expected relatedness
#' based on pedigree and (2) realized relatedness based on marker data. See
#' 'Details'.  The function uses as first argument an object of class
#' \code{gpData}. An argument \code{ret} controls the type of relatedness
#' coefficient.
#'
#'
#' \bold{Pedigree based relatedness (return arguments \code{"add"},
#' \code{"kin"}, \code{"dom"}, and \code{"gam"})}
#'
#' Function \code{kin} provides different types of measures for pedigree based
#' relatedness. An element \code{pedigree} must be available in the object of
#' class \code{gpData}. In all cases, the first step is to build the gametic
#' relationship. The gametic relationship is of order 2\eqn{n} as each
#' individual has two alleles (e.g. individual \eqn{A} has alleles \eqn{A1} and
#' \eqn{A2}). The gametic relationship is defined as the matrix of
#' probabilities that two alleles are identical by descent (IBD).  Note that
#' the diagonal elements of the gametic relationship matrix are 1. The
#' off-diagonals of individuals with unknown or unrelated parents in the
#' pedigree are 0. If \code{ret="gam"} is specified, the gametic relationship
#' matrix constructed by pedigree is returned.
#'
#' The gametic relationship matrix can be used to construct other types of
#' relationship matrices. If \code{ret="add"}, the additive numerator
#' relationship matrix is returned. The additive relationship of individuals A
#' (alleles \eqn{A1,A2}) and B (alleles \eqn{B1,B2}) is given by the entries of
#' the gametic relationship matrix \deqn{0.5\cdot \left[(A1,B1) + (A1,B2) +
#' (A2,B1) + (A2,B2)\right],}{0.5*[(A1,B1) + (A1,B2) + (A2,B1) + (A2,B2)],}
#' where \eqn{(A1,B1)} denotes the element [A1,B1] in the gametic relationship
#' matrix. If \code{ret="kin"}, the kinship matrix is returned which is half of
#' the additive relationship matrix.
#'
#' If \code{ret="dom"}, the dominance relationship matrix is returned. The
#' dominance relationship matrix between individuals A (\eqn{A1,A2}) and B
#' (\eqn{B1,B2}) in case of no inbreeding is given by \deqn{\left[(A1,B1) \cdot
#' (A2,B2) + (A1,B2) \cdot (A2,B1)\right],}{[(A1,B1) * (A2,B2) + (A1,B2) *
#' (A2,B1)],} where \eqn{(A1,C1)} denotes the element [A1,C1] in the gametic
#' relationship matrix.
#'
#' \bold{Marker based relatedness (return arguments
#' \code{"realized"},\code{"realizedAB"}, \code{"sm"}, and \code{"sm-smin"})}
#'
#' Function \code{kin} provides different types of measures for marker based
#' relatedness. An element \code{geno} must be available in the object of class
#' \code{gpData}. Furthermore, genotypes must be coded by the number of copies
#' of the minor allele, i.e. function \code{codeGeno} must be applied in
#' advance.
#'
#' If \code{ret="realized"}, the realized relatedness between individuals is
#' computed according to the formulas in Habier et al. (2007) or vanRaden
#' (2008) \deqn{U = \frac{ZZ'}{2\sum p_i(1-p_i)}}{ZZ'/(2\sum pi(1-pi))} where
#' \eqn{Z=W-P}, \eqn{W} is the marker matrix, \eqn{P} contains the allele
#' frequencies multiplied by 2, \eqn{p_i}{pi} is the allele frequency of marker
#' \eqn{i}, and the sum is over all loci.
#'
#' If \code{ret="realizedAB"}, the realized relatedness between individuals is
#' computed according to the formula in Astle and Balding (2009) \deqn{U =
#' \frac{1}{M} \sum \frac{(w_i-2p_i)(w_i-2p_i)'}{2p_i(1-p_i)}}{1/M
#' sum((wi-2pi)(wi-2pi)'/(2pi(1-pi)))} where \eqn{w_i}{wi} is the marker
#' genotype, \eqn{p_i}{pi} is the allele frequency at marker locus \eqn{i}, and
#' \eqn{M} is the number of marker loci, and the sum is over all loci.
#'
#' If \code{ret="sm"}, the realized relatedness between individuals is computed
#' according to the simple matching coefficient (Reif et al. 2005). The simple
#' matching coefficient counts the number of shared alleles across loci. It can
#' only be applied to homozygous inbred lines, i.e. only genotypes 0 and 2. To
#' account for loci that are alike in state but not identical by descent (IBD),
#' Hayes and Goddard (2008) correct the simple matching coefficient by the
#' minimum of observed simple matching coefficients
#' \deqn{\frac{s-s_{min}}{1-s_{min}}}{s-smin/(1-smin)} where \eqn{s} is the
#' matrix of simple matching coefficients. This formula is used with argument
#' \code{ret="sm-smin"}.
#'
#' If \code{ret="gaussian"}, the euklidian distances \code{distEuk} for all
#' individuals are calculated. The values of \code{distEuk} are than used to
#' calculate similarity coefficients between the individuals with
#' \code{exp(distEuk^2/numMarker)}. Be aware that this relationship matrix
#' scales theoretically between 0 and 1!
#'
#' @param gpData object of class \code{gpData}
#' @param ret \code{character}. The type of relationship matrix to be returned.
#' See 'Details'.
#' @param DH \code{logical} vector of length \eqn{n}. \code{TRUE} or 1 if
#' individual is a doubled-haploid (DH) line and \code{FALSE} or 0 otherwise.
#' This option is only used, if \code{ret} argument is \code{"add"} or
#' \code{"kin"}.
#' @param maf \code{numeric} vector of length equal the number of markers.
#' Supply values for the \eqn{p_i}{pi} of each marker, which were used to
#' correct the allele counts in \code{ret="realized"} and
#' \code{ret="realizedAB"}. If not specified, \eqn{p_i}{pi} equals the minor
#' allele frequency of each locus.
#' @param selfing \code{numeric} vector of length \eqn{n}. It is used as the
#' number of selfings of an recombinant inbred line individual. Be awere, that
#' this should only be used for single seed descendants This option is only
#' used, if \code{ret} argument is \code{"add"} or \code{"kin"}.
#' @param lambda \code{numeric} bandwidth parameter for the gaussian kernel.
#' Only used for calculating the gaussian kernel.
#' @param P \code{numeric} matrix of the same dimension as \code{geno} of the
#' \code{gpData} object. This option can be used for own allelefrequencies of
#' different groups in the genotypes.
#' @param cores \code{numeric}. Here you can specify the number of cores you
#' like to use.
#' @return An object of class "relationshipMatrix".
#' @author Valentin Wimmer and Theresa Albrecht, with contributions by Yvonne
#' Badke
#' @seealso \code{\link{plot.relationshipMatrix}}
#' @references
#'
#' Habier D, Fernando R, Dekkers J (2007). The Impact of Genetic Relationship
#' information on Genome-Assisted Breeding Values. Genetics, 177, 2389 -- 2397.
#'
#' vanRaden, P. (2008). Efficient methods to compute genomic predictions.
#' Journal of Dairy Science, 91:4414 -- 4423.
#'
#' Astle, W., and D.J. Balding (2009). Population Structure and Cryptic
#' Relatedness in Genetic Association Studies. Statistical Science, 24(4), 451
#' -- 471.
#'
#' Reif, J.C.; Melchinger, A. E. and Frisch, M. Genetical and mathematical
#' properties of similarity and dissimilarity coefficients applied in plant
#' breeding and seed bank management. Crop Science, January-February 2005, vol.
#' 45, no. 1, p. 1-7.
#'
#' Rogers, J., 1972 Measures of genetic similarity and genetic distance. In
#' Studies in genetics VII, volume 7213. Univ. of Texas, Austin
#'
#' Hayes, B. J., and M. E. Goddard. 2008. Technical note: Prediction of
#' breeding values using marker derived relationship matrices. J. Anim. Sci. 86
#' @examples
#'
#'
#' # =========================
#' # (1) pedigree based relatedness
#' # =========================
#' \dontrun{
#' library(synbreedData)
#' data(maize)
#' K <- kin(maize, ret = "kin")
#' plot(K)
#' }
#'
#' # =========================
#' # (2) marker based relatedness
#' # =========================
#' \dontrun{
#' data(maize)
#' U <- kin(codeGeno(maize), ret = "realized")
#' plot(U)
#' }
#'
#'
#' ### Example for Legarra et al. (2009), J. Dairy Sci. 92: p. 4660
#' id <- 1:17
#' par1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 5, 7, 9, 11, 4, 13, 13)
#' par2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 6, 8, 10, 12, 11, 15, 14)
#' ped <- create.pedigree(id, par1, par2)
#' gp <- create.gpData(pedigree = ped)
#'
#' # additive relationship
#' A <- kin(gp, ret = "add")
#' # dominance relationship
#' D <- kin(gp, ret = "dom")
#' @export kin
#' @importFrom stats dist
#' @importFrom utils sessionInfo
#' @importFrom graphics plot
#'
kin <- function(gpData, ret = c("add", "kin", "dom", "gam", "realized", "realizedAB", "sm", "sm-smin", "gaussian"),
                DH = NULL, maf = NULL, selfing = NULL, lambda = 1, P = NULL, cores = 1) {
  ret <- match.arg(ret, choices = c("add", "kin", "dom", "gam", "realized", "realizedAB", "sm", "sm-smin", "gaussian"), several.ok = FALSE)
  NAs <- FALSE
  if (ret %in% c("realized", "realizedAB", "sm", "sm-smin")) if (any(is.na(gpData$geno))) NAs <- TRUE
  # (1) expected relatedness
  if (ret %in% c("add", "kin", "dom", "gam")) {
    # check for 'gpData'
    if (any(class(gpData) == "gpData")) {
      if (is.null(gpData$pedigree)) {
        stop("no pedigree found")
      } else {
        ped <- gpData$pedigree
      }
    } else if (any(class(gpData) == "pedigree")) {
      ped <- gpData
    }

    # number of ids
    n <- nrow(ped)
    if (is.null(DH)) DH <- rep(0, n)
    if (!is.null(DH) & (length(DH) != n)) stop("DH must have same length as pedigree")
    if (!is.null(selfing) & (length(selfing) != n)) stop("DH must have same length as pedigree")

    if (ret %in% c("add", "kin")) {
      IDplus <- unique(c(0, ped$Par1[!ped$Par1 %in% ped$ID], ped$Par2[!ped$Par2 %in% ped$ID]))
      if (length(IDplus) > 1) warning("There are parents in the pedigree, which are not coded as ancestors in the ID column!")
      A <- matrix(data = 0, nrow = n + length(IDplus), ncol = n + length(IDplus))
      if (is.null(selfing)) selfing <- rep(0, ncol(A)) else selfing <- c(rep(0, length(IDplus)), selfing)
      if (is.numeric(c(IDplus, ped$ID))) stop("You use only numbers as identifier of your individuals! This causes trouble!")
      names(selfing) <- colnames(A) <- rownames(A) <- c(IDplus, ped$ID)
      A[IDplus, IDplus] <- diag(length(IDplus))
      A[1, 1] <- 0
      for (i in ped$ID) {
        cnt <- match(i, ped$ID)
        cnt2 <- cnt + length(IDplus)
        A[i, 1:cnt2] <- A[1:cnt2, i] <- (A[ped[cnt, "Par1"], 1:cnt2] + A[ped[cnt, "Par2"], 1:cnt2]) * .5
        if (DH[cnt] != 0) {
          A[i, i] <- 2
        } else {
          A[i, i] <- 2 - .5^selfing[i] * (1 - .5 * A[ped[cnt, "Par1"], ped[cnt, "Par2"]])
        } # Schoenleben et al., unpublished
      }
      A <- A[-c(1:length(IDplus)), -c(1:length(IDplus))]
    } else if (ret %in% c("gam", "dom")) {
      # set up extended pedigree
      ID <- rep(seq_along(ped$ID), each = 2)
      par1 <- pmatch(ped$Par1, ped$ID, nomatch = 0, duplicates.ok = TRUE)
      par2 <- pmatch(ped$Par2, ped$ID, nomatch = 0, duplicates.ok = TRUE)

      # set up gametic pedigree data.frame
      gamMat <- matrix(data = 0, nrow = n * 2, ncol = 3, byrow = FALSE)
      gamMat[, 1] <- ID
      # loop over ID
      for (i in 1:n) {
        par1gam <- par1[i]
        par2gam <- par2[i]
        j <- (i - 1) * 2 + 1
        k <- j + 1
        #  parents of male genome contribution
        if (par1gam > 0) {
          gamMat[j, 2] <- (par1gam - 1) * 2 + 1
          gamMat[j, 3] <- (par1gam - 1) * 2 + 2
        }
        #  parents of female genome contribution
        if (par2gam > 0) {
          gamMat[k, 2] <- (par2gam - 1) * 2 + 1
          gamMat[k, 3] <- (par2gam - 1) * 2 + 2
        }
      } # end of loop over ID

      #  Build Gametic Relationship
      ngam <- 2 * n
      DHgam <- rep(DH, each = 2)
      G <- diag(ngam)
      dimnames(G) <- list(paste(rep(ped$ID, each = 2), rep(1:2, times = n), sep = "_"), paste(rep(ped$ID, each = 2), rep(1:2, times = n), sep = "_"))
      # set inbreed coefficients of DHs on 1
      G[cbind((1:ngam) * DHgam, ((1:ngam) + c(1, -1)) * DHgam)] <- 1

      # caluclate gametic relationship
      # loop over gamets
      for (i in 1:(ngam - 1 - DHgam[2 * n])) {
        ip <- i + 1 + (DHgam * rep(c(1, 0), ngam))[i]
        for (j in ip:ngam) {
          if (gamMat[j, 2] > 0) {
            x <- 0.5 * (G[i, gamMat[j, 2]] + G[i, gamMat[j, 3]])
            G[i, j] <- G[j, i] <- x
          }
        }
      } # end of loop over gamets

      # calculate dominance relationship
      if (ret == "dom") {
        D <- matrix(data = NA, nrow = n, ncol = n)
        dimnames(D) <- list(ped$ID, ped$ID)

        # set up D matrix
        # loop over individuals
        for (i in 1:n) {
          ka <- (i - 1) * 2 + 1
          for (j in i:n) {
            kb <- (j - 1) * 2 + 1
            dab <- (G[ka, kb] * G[ka + 1, kb + 1] + G[ka + 1, kb] * G[ka, kb + 1]) #* (1-G[ka,ka+1])*(1-G[kb,kb+1])
            # acoount for inbreeding
            # dominance = 0 if Fi=1
            D[i, j] <- D[j, i] <- dab
          }
        } # end of loop over individuals
      } # end of if
    }

    # set return matrices
    if (ret == "add") kmat <- A
    if (ret == "dom") kmat <- D
    if (ret == "kin") kmat <- A / 2
    if (ret == "gam") kmat <- G

    attr(kmat, "SNPs") <- NULL
  }

  # (2) realized relatedness

  if (ret == "realized") { # former method vanRaden
    # extract information from arguments
    if (any(class(gpData) == "gpData")) {
      if (!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'")
      W <- gpData$geno
      if (!is.null(maf) & length(maf) != ncol(gpData$geno)) stop("minor allele frequency not provided for all markers")
    } else {
      stop("object is not of class 'gpData'")
    }

    # W supposed to be coded with 0,1,2
    n <- nrow(W)
    p <- ncol(W)

    # use user-supplied values for maf
    # or, otherwise 2* minor allele frequency as expectation
    if (is.null(maf)) {
      maf <- colMeans(W, na.rm = TRUE)
    }

    if (is.null(P)) {
      p <- matrix(rep(maf, each = n), ncol = p)
    } else if (!all.equal(dim(P), dim(W))) {
      stop("wrong dimension of the matrix P")
    } else {
      p <- as.matrix(P)
    }
    # compute realized relationship matrix G
    Z <- W - p
    if (NAs) Z[is.na(Z)] <- 0
    U <- tcrossprod(Z)
    U <- 2 * U / (sum(maf * (2 - maf)))

    kmat <- U
    if (!is.null(P)) attr(kmat, "P") <- P
    attr(kmat, "alleleFrequencies") <- maf
    attr(kmat, "expectedMAX") <- 2 * sum((2 - maf)**2) / (sum(maf * (2 - maf)))
    attr(kmat, "SNPs") <- colnames(gpData$geno)
  }

  if (ret == "realizedAB") { # based an Astle & Balding (2009)

    # extract information from arguments
    if (any(class(gpData) == "gpData")) {
      if (!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'")
      W <- gpData$geno
      if (!is.null(maf) & length(maf) != ncol(gpData$geno)) stop("minor allele frequency not provided for all markers")
    }
    else {
      stop("object is not of class 'gpData'")
    }

    # W supposed to be coded with 0,1,2
    n <- nrow(W)
    p <- ncol(W)

    # use user-supplied values for maf
    # or, otherwise 2* minor allele frequency as expectation
    if (is.null(maf)) {
      maf <- colMeans(W, na.rm = TRUE)
    }

    pq2 <- .5 * maf * (2 - maf)
    # compute realized relationship matrix U
    W <- sweep(W, 2, maf)
    W <- sweep(W, 2, sqrt(pq2), "/")
    if (NAs) W[is.na(W)] <- 0
    U <- tcrossprod(W) / p
    kmat <- U
    attr(kmat, "alleleFrequencies") <- maf
    attr(kmat, "markerVariances") <- pq2
    attr(kmat, "SNPs") <- colnames(gpData$geno)
  }

  if (ret %in% c("sm", "sm-smin")) { # simple matchin coefficient (only for homozygous inbreed lines)

    # extract information from arguments
    if (any(class(gpData) == "gpData")) {
      if (!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'")
      marker <- gpData$geno
    }
    else {
      stop("object is not of class 'gpData'")
    }

    # code marker to -1/0/1 from 0,1,2
    marker <- marker - (max(marker, na.rm = TRUE) - 1)
    m <- ncol(marker)
    if (NAs) marker[is.na(marker)] <- 0
    s <- (tcrossprod(marker) + m) / (2 * m)

    if (ret == "sm-smin") {
      smin <- min(s, na.rm = TRUE)
      s <- (s - smin) / (1 - smin)
      attr(s, "min") <- smin
    }
    kmat <- 2 * s
    attr(kmat, "SNPs") <- colnames(gpData$geno)
  }

  if (ret == "gaussian") { # Euclidian distance with Gaussian
    if (any(class(gpData) == "gpData")) {
      if (!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'")
      marker <- gpData$geno
    } else {
      stop("object is not of class 'gpData'")
    }

    marker <- scale(marker, center = TRUE, scale = TRUE)
    Dist <- (as.matrix(dist(marker, method = "euclidean"))**2) / ncol(marker)
    kmat <- exp(-lambda * Dist)

    attr(kmat, "SNPs") <- colnames(gpData$geno)
  }
  attr(kmat, "info") <- paste("This relationshipMatrix was calculated by synbreed version", sessionInfo()$otherPkgs$synbreed$Version)
  attr(kmat, "type") <- ret
  class(kmat) <- c("relationshipMatrix", "matrix")
  return(kmat)
}

#' Simulation of a field trial with single trait
#'
#' Simulates observations from a field trial using an animal model. The field
#' trial consists of multiple locations and randomized complete block design
#' within locations.  A single quantitative trait is simulated according to the
#' model \code{Trait ~ id(A) + block + loc + e}.
#'
#' Either \code{pedigree} or \code{A} must be specified. If \code{pedigree} is
#' given, pedigree information is used to set up numerator relationship matrix
#' with function \code{kinship}. If unrelated individuals should be used for
#' simulation, use identity matrix for \code{A}. True breeding values for
#' \eqn{N} individuals is simulated according to following distribution
#' \deqn{tbv \sim N(0,\bf A \sigma_a^2)}{tbv = chol(A)*sigma2a*rnorm(N,0,1)}
#' Observations are simulated according to \deqn{y \sim N(mu + tbv + block
#' +loc,\sigma^2_e)}{} If no location or block effects should appear, use
#' \code{sigma2l=0} and/or \code{sigma2b=0}.
#'
#' @param pedigree object of class "pedigree"
#' @param A object of class "relationshipMatrix"
#' @param mu \code{numeric}; Overall mean of the trait.
#' @param vc list containing the variance components. \code{vc} consists of
#' elements \code{sigma2e}, \code{sigma2a}, \code{sigma2l}, \code{sigma2b} with
#' the variance components of the residual, the additive genetic effect, the
#' location effect and the block effect.
#' @param Nloc \code{numeric}. Number of locations in the field trial.
#' @param Nrepl \code{Numeric}. Number of complete blocks within location.
#' @return A \code{data.frame} with containing the simulated values for trait
#' and the following variables \item{ID}{Factor identifying the individuals.
#' Names are extracted from \code{pedigree} or rownames of \code{A}}
#' \item{Loc}{Factor for Location} \item{Block}{Factor for Block within
#' Location} \item{Trait}{Trait observations} \item{TBV}{Simulated values for
#' true breeding values of individuals} Results are sorted for locations within
#' individuals.
#' @author Valentin Wimmer
#' @seealso \code{\link{simul.pedigree}}
#' @examples
#'
#' \dontrun{
#' ped <- simul.pedigree(gener = 5)
#' varcom <- list(sigma2e = 25, sigma2a = 36, sigma2l = 9, sigma2b = 4)
#' # field trial with 3 locations and 2 blocks within locations
#' data.simul <- simul.phenotype(ped, mu = 10, vc = varcom, Nloc = 3, Nrepl = 2)
#' head(data.simul)
#' # analysis of variance
#' anova(lm(Trait ~ ID + Loc + Loc:Block, data = data.simul))
#' }
#'
#' @export simul.phenotype
#' @importFrom stats lm
#' @importFrom stats rnorm step
#'
simul.phenotype <- function(pedigree = NULL, A = NULL, mu = 100, vc = NULL, Nloc = 1, Nrepl = 1) {
  if (is.null(A) & is.null(pedigree)) step("either 'pedigree' or 'A' must be given")
  # create A matrix if missing
  if (is.null(A)) A <- kin(pedigree, ret = "add")
  if (is.null(vc)) stop("missing variance components")

  # read data out of arguments
  N <- nrow(A)
  n <- N * Nloc * Nrepl
  if (!is.null(vc$sigma2e)) {
    sigmae <- sqrt(vc$sigma2e)
  } else {
    sigmae <- 1
  }
  if (!is.null(vc$sigma2a)) {
    sigmaa <- sqrt(vc$sigma2a)
  } else {
    sigmaa <- 0
  }
  if (!is.null(vc$sigma2l)) {
    sigmal <- sqrt(vc$sigma2l)
  } else {
    sigmal <- 0
  }
  if (!is.null(vc$sigma2b)) {
    sigmab <- sqrt(vc$sigma2b)
  } else {
    sigmab <- 0
  }

  namesA <- rownames(A)
  if (is.null(rownames(A))) namesA <- 1:N

  # initialize data

  ID <- rep(namesA, each = Nrepl * Nloc)
  Loc <- rep(1:Nloc, length.out = n, each = Nrepl)
  Block <- rep(1:Nrepl, lengt.out = n)

  # as matrix
  A <- matrix(A, nrow = N)


  # simulate data for contribution parts of phenotypic values

  # true breeding values
  # tbv <- rmvnorm(1,rep(0,N),A*sigmaa^2)
  tbv <- sigmaa * (chol(A) %*% rnorm(N, 0, 1))
  tbv <- rep(tbv, each = Nloc * Nrepl)
  # location effect
  locEff <- rnorm(Nloc, 0, sigmal)
  locEff <- rep(locEff, length.out = n, each = Nrepl)
  # block effects
  blockEff <- rnorm(Nrepl, 0, sigmab)
  blockEff <- rep(blockEff, length.out = n)
  # residual
  residual <- rnorm(n, 0, sigmae)

  # simulate phenotypic value
  Trait <- mu + tbv + locEff + blockEff + residual

  # combine to a data.frame
  ret <- data.frame(ID = factor(ID), Loc = factor(Loc), Block = factor(Block), Trait = Trait, TBV = tbv)
  return(ret)
}

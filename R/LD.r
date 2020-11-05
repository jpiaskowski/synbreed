#' Pairwise LD between markers
#'
#' Estimate pairwise Linkage Disequilibrium (LD) between markers measured as
#' \eqn{r^2}{r2} using an object of class \code{gpData}. For the general case,
#' a gateway to the software PLINK (Purcell et al. 2007) is established to
#' estimate the LD. A within-R solution is only available for marker data with
#' only 2 genotypes, i.e. homozgous inbred lines. Return value is an object of
#' class \code{LDdf} which is a \code{data.frame} with one row per marker pair
#' or an object of class \code{LDMat} which is a \code{matrix} with all marker
#' pairs. Additionally, the euclidian distance between position of markers is
#' computed and returned.
#'
#' The function \code{write.plink} is called to prepare the input files and the
#' script for PLINK. The executive PLINK file \code{plink.exe} must be
#' available (e.g. in the working directory or through path variables). The
#' function \code{pairwiseLD} calls PLINK and reads the results. The evaluation
#' is performed separately for every chromosome. The measure for LD is
#' \eqn{r^2}{r2}. This is defined as \deqn{D= p_{AB} - p_Ap_B }{D = p(AB) -
#' p(A)*p(B)} and \deqn{r^2=\frac{D^2}{p_Ap_Bp_ap_b}}{r2=D2/p(A)p(B)p(a)p(b)}
#' where \eqn{p_{AB}}{p(AB)} is defined as the observed frequency of haplotype
#' \eqn{AB}, \eqn{p_A=1-p_a} and \eqn{p_B=1-p_b} the observed frequencies of
#' alleles \eqn{A} and \eqn{B}. If the number of markers is high, a threshold
#' for the LD can be used to thin the output. In this case, only pairwise LD
#' above the threshold is reported (argument \code{--ld-window-r2 in PLINK}).
#'
#' Default PLINK options used --no-parents --no-sex --no-pheno --allow-no-sex
#' --ld-window p --ld-window-kb 99999
#'
#' @param gpData object of class \code{gpData} with elements \code{geno} and
#' \code{map}
#' @param chr \code{numeric} scalar or vector. Return value is a list with
#' pairwise LD of all markers for each chromosome in \code{chr}.
#' @param type \code{character}. Specifies the type of return value (see
#' 'Value').
#' @param use.plink \code{logical}. Should the software PLINK be used for the
#' computation?
#' @param ld.threshold \code{numeric}. Threshold for the LD to thin the output.
#' Only pairwise LD>\code{ld.threshold} is reported when PLINK is used. This
#' argument can only be used for \code{type="data.frame"}.
#' @param ld.window \code{numeric}. Window size for pairwise differences which
#' will be reported by PLINK (only for \code{use.plink=TRUE}; argument
#' \code{--ld-window-kb} in PLINK) to thin the output dimensions. Only SNP
#' pairs with a distance < \code{ld.window} are reported (default = 99999).
#' @param rm.unmapped \code{logical}. Remove markers with unknown postion in
#' \code{map} before using PLINK?
#' @param cores \code{numeric}. Here you can specify the number of cores you
#' like to use.
#' @return For \code{type="data.frame"} an object of class \code{LDdf} with one
#' element for each chromosome is returned. Each element is a \code{data.frame}
#' with columns \code{marker1}, \code{marker2}, \code{r2} and \code{distance}
#' for all \eqn{p(p-1)/2} marker pairs (or thinned, see 'Details').
#'
#' For \code{type="matrix"} an object of class \code{LDmat} with one element
#' for each chromosome is returned. Each element is a list of 2: a \eqn{p
#' \times p}{p x p} \code{matrix} with pairwise LD and the corresponding \eqn{p
#' \times p}{p x p} \code{matrix} with pairwise distances.
#' @author Valentin Wimmer
#' @seealso \code{\link{LDDist}}, \code{\link{LDMap}}
#' @references Hill WG, Robertson A (1968). Linkage Disequilibrium in Finite
#' Populations. Theoretical and Applied Genetics, 6(38), 226 - 231.
#'
#' Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller
#' J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) PLINK: a toolset for
#' whole-genome association and population-based linkage analysis. American
#' Journal of Human Genetics, 81.
#' @examples
#'
#' \dontrun{
#' library(synbreedData)
#' data(maize)
#' maizeC <- codeGeno(maize)
#' maizeLD <- pairwiseLD(maizeC, chr = 1, type = "data.frame")
#' }
#'
#' @export pairwiseLD
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores makeCluster parLapply stopCluster mclapply
#' @importFrom stats cor dist
#' @importFrom utils read.table
#'
pairwiseLD <- function(gpData, chr = NULL, type = c("data.frame", "matrix"), use.plink = FALSE,
                       ld.threshold = 0, ld.window = 99999, rm.unmapped = TRUE, cores = 1) {
  multiLapply <- function(x, y, ..., mc.cores = 1) {
    if (.Platform$OS.type == "windows" & cores > 1) {
      cl <- makeCluster(min(mc.cores, detectCores()))
      registerDoParallel(cl)
      parLapply(cl, x, y, ...)
      stopCluster(cl)
    } else {
      mclapply(x, y, ..., mc.cores = mc.cores)
    }
  }
  multiCor <- function(x, use = "everything", method = c("pearson", "kendall", "spearman"), cores = 1) {
    if (cores == 1) {
      cor(x, use = use, method = method)
    } else {
      method <- match.arg(method)
      ncolX <- ncol(x)
      namesX <- colnames(x)
      mat <- matrix(NA, nrow = ncolX, ncol = ncolX)
      nBlocks <- 1:50
      nBlocks <- which.min((nBlocks * (nBlocks + 1) * .5) %% cores)
      if (nBlocks == 1) {
        nBlocks <- 1:50
        nBlocks <- which.min((nBlocks * (nBlocks + 1) * .5) %% (cores - 1))
      }
      blocks <- unique(t(apply(cbind(rep(1:nBlocks, nBlocks), rep(1:nBlocks, each = nBlocks)), 1, sort)))
      groups <- split(1:ncolX, sort(rep(1:nBlocks, ceiling(ncolX / nBlocks))[1:ncolX]))
      cl <- makeCluster(min(cores, detectCores()))
      registerDoParallel(cl)
      cors <- list()
      res <- foreach(i = 1:nrow(blocks)) %dopar% {
        cors[[i]] <- cor(x[, groups[[blocks[i, 1]]]], x[, groups[[blocks[i, 2]]]], use = use, method = method)
      }
      stopCluster(cl)
      for (i in 1:nrow(blocks)) {
        mat[groups[[blocks[i, 1]]], groups[[blocks[i, 2]]]] <- res[[i]]
        mat[groups[[blocks[i, 2]]], groups[[blocks[i, 1]]]] <- t(res[[i]])
      }
      return(mat)
    }
  }

  # catch errors
  if (is.null(gpData$geno)) stop("no genotypic data available")
  if (!gpData$info$codeGeno) stop("use function 'codeGeno' before")
  if (is.null(gpData$map)) stop("no map information available")
  type <- match.arg(type)
  if (ld.threshold > 0 & type == "matrix") warning("'ld.threshold not used for type='matrix'")
  if (any(is.na(gpData$geno))) stop("no missing values allowed, try to impute using 'codeGeno'")

  # extract information from gpData  (only use mapped markers)
  if (rm.unmapped) {
    mapped <- !(is.na(gpData$map$chr) | is.na(gpData$map$pos))
    gpData$geno <- gpData$geno[, mapped]
    gpData$map <- gpData$map[mapped, ]
  } else {
    gpData$map$chr <- as.character(gpData$map$chr)
    gpData$map$chr[is.na(gpData$map$pos)] <- "NA"
    gpData$map$chr <- as.factor(gpData$map$chr)
  }
  linkageGroup <- as.character(gpData$map$chr)
  pos <- gpData$map$pos
  names(pos) <- rownames(gpData$map)

  # select chromosomes if 'chr' is specified
  lg <- unique(linkageGroup)
  if (!is.null(chr)) {
    lg <- chr
    if (any(chr == "all")) stop("option chr='all' not yet possible") # linkageGroup <- rep("all",length(linkageGroup))
    # NOTE: positions must be ordered consequtively within chromosomes
  }

  # initialize return LD value data.frame list
  retList <- list()
  # initialize return LD value matrix list
  retMat <- list()

  # loop over all chromosomes (linkage groups)
  for (i in 1:length(lg)) {
    if (use.plink) { # i.e. if there are 3 genotypes
      # call PLINK to compute the LD as r2
      sel <- rownames(gpData$map)[gpData$map$chr != lg[i]]
      gpTEMP <- discard.markers(gpData, which = sel)
      pre <- paste("chr", lg[i], sep = "")
      write.plink(gpTEMP, type = type, ld.threshold = ld.threshold, ld.window = ld.window, prefix = pre)
      system(paste("plink --script ", pre, "plinkScript.txt", sep = ""))
      # distances between markers
      if (type == "matrix") {
        distance <- as.matrix(dist(pos[linkageGroup == lg[i]], diag = FALSE, upper = FALSE))
      }
      # read data from PLINK
      if (type == "matrix") {
        ld.r2 <- as.matrix(read.table(paste(pre, ".ld", sep = "")))
        colnames(distance) <- rownames(distance) <- colnames(ld.r2) <- rownames(ld.r2) <- names(pos)[linkageGroup == lg[i]]
        ld.r <- sqrt(ld.r2)
      }
      if (type == "data.frame") {
        ld.r2.df.plink <- read.table(paste(pre, ".ld", sep = ""), header = TRUE, stringsAsFactors = FALSE)
        # distance <- abs(pos[ld.r2.df.plink$SNP_A]-pos[ld.r2.df.plink$SNP_B])
        ld.r2.df <- with(ld.r2.df.plink, data.frame(marker1 = SNP_A, marker2 = SNP_B, r2 = R2, dist = abs(BP_A - BP_B), stringsAsFactors = FALSE))
      }
    } # end if(use.plink)
    else { # i.e. if there are 2 genotypes (e.g. DH lines)
      # read information from data
      markeri <- gpData$geno[, linkageGroup == lg[i]]
      p <- ncol(markeri)
      mn <- colnames(markeri)
      posi <- pos[linkageGroup == lg[i]]
      ld.r <- multiCor(markeri, method = "spearman", use = "pairwise.complete.obs", cores = cores)
      ld.r2 <- ld.r**2
      if (type == "data.frame") {
        ld.ri <- ld.r[lower.tri(ld.r)]
        # index vectors for LD data.frame
        rowi <- rep(1:p, times = (p:1) - 1)
        coli <- p + 1 - sequence(1:(p - 1))
        coli <- coli[length(coli):1]
        # distance between markers
        disti <- abs(posi[rowi] - posi[coli])
        ld.r2.df <- data.frame(
          marker1 = mn[rowi],
          marker2 = mn[coli],
          r = ld.ri,
          r2 = ld.ri**2,
          dist = disti,
          stringsAsFactors = FALSE
        )
        if (lg[i] == "NA") {
          ld.rini <- multiCor(gpData$geno[, linkageGroup != lg[i]], markeri, method = "spearman", use = "pairwise.complete.obs", cores = cores)
          ld.r2.dfini <- data.frame(
            marker1 = rep(colnames(ld.rini), each = nrow(ld.rini)),
            marker2 = rep(rownames(ld.rini), ncol(ld.rini)),
            r = as.numeric(ld.rini),
            r2 = as.numeric(ld.rini**2),
            dist = rep(NA, ncol(ld.rini) * nrow(ld.rini)),
            stringsAsFactors = FALSE
          )
          ld.r2.df <- rbind(ld.r2.df, ld.r2.dfini)
        }
      }
      if (type == "matrix") {
        # matrix of distances
        distance <- as.matrix(dist(pos[linkageGroup == lg[i]], diag = FALSE, upper = FALSE))
        colnames(distance) <- rownames(distance) <- colnames(gpData$geno)[linkageGroup == lg[i]]
      }
    }
    # create dataset with information from above in a data.frame
    if (type == "data.frame") retList[[i]] <- ld.r2.df
    # and as a matrix
    if (type == "matrix") retMat$LD[[i]] <- ld.r2 # omit lower/upper triangle?
    if (type == "matrix") retMat$distance[[i]] <- distance
    if (type == "matrix") retMat$LDcor[[i]] <- ld.r # omit lower/upper triangle?
  }
  if (type == "data.frame") names(retList) <- paste("chr", lg, sep = "_")
  if (type == "matrix") names(retMat$LD) <- names(retMat$distance) <- paste("chr", lg, sep = "_")

  # return values
  # if(type=="both") return(list(dataFrame=retList,matrix=retMat))
  if (type == "data.frame") {
    class(retList) <- "LDdf"
    return(retList)
  }
  if (type == "matrix") {
    class(retMat) <- "LDmat"
    return(retMat)
  }
}

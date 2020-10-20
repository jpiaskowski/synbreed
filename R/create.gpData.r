# read genomic prediction data


#' Create genomic prediction data object
#'
#' This function combines all raw data sources in a single, unified data object
#' of class \code{gpData}. This is a \code{list} with elements for phenotypic,
#' genotypic, marker map, pedigree and further covariate data. All elements are
#' optional.
#'
#' The class \code{gpData} is designed to provide a unified framework for data
#' related to genomic prediction analysis. Every data source can be omitted. In
#' this case, the corresponding argument must be \code{NULL}. By default
#' (argument \code{reorderMap}), markers in \code{geno} are ordered by their
#' position in \code{map}. Individuals are ordered in alphabetical order.
#'
#' An object of class \code{gpData} can contain different subsets of
#' individuals or markers in the elements \code{pheno}, \code{geno} and
#' \code{pedigree}. In this case the \code{id} in \code{covar} comprises all
#' individuals that either appear in \code{pheno}, \code{geno} and
#' \code{pedigree}. Two additional columns in \code{covar} named
#' \code{phenotyped} and \code{genotyped} are automatically generated to
#' identify individuals that appear in the corresponding \code{gpData} object.
#'
#' @param pheno \code{data.frame} with individuals organized in rows and traits
#' organized in columns. For unrepeated measures unique \code{rownames} should
#' identify individuals. For repeated measures, the first column identifies
#' individuals and a second column indicates repetitions (see also argument
#' \code{repeated}).
#' @param geno \code{matrix} with individuals organized in rows and markers
#' organized in columns. Genotypes could be coded arbitrarily. Missing values
#' should be coded as \code{NA}. Colums or rows with only missing values not
#' allowed. Unique \code{rownames} identify individuals and unique
#' \code{colnames} markers. If no \code{rownames} are available, they are taken
#' from element \code{pheno} (if available and if dimension matches).  If no
#' \code{colnames} are used, the \code{rownames} of \code{map} are used if
#' dimension matches.
#' @param map \code{data.frame} with one row for each marker and two columns
#' (named \code{chr} and \code{pos}). First columns gives the chromosome
#' (\code{numeric} or \code{character} but not \code{factor}) and second column
#' the position on the chromosome in centimorgan or the physical distance
#' relative to the reference sequence in basepairs. Unique \code{rownames}
#' indicate the marker names which should match with marker names in
#' \code{geno}. Note that order and number of markers must not be identical
#' with the order in \code{geno}. If this is the case, gaps in the map are
#' filled with \code{NA} to ensure the same number and order as in element
#' \code{geno} of the resulting \code{gpData} object.
#' @param pedigree Object of class \code{pedigree}.
#' @param family \code{data.frame} assigning individuals to families with names
#' of individuals in \code{rownames} This information could be used for
#' replacing of missing values with function \code{codeGeno}.
#' @param covar \code{data.frame} with further covariates for all individuals
#' that either appear in \code{pheno}, \code{geno} or \code{pedigree$ID}, e.g.
#' sex or age. \code{rownames} must be specified to identify individuals.
#' Typically this element is not specified by the user.
#' @param reorderMap \code{logical}. Should markers in \code{geno} and
#' \code{map} be reordered by chromosome number and position within chromosome
#' according to \code{map} (default = \code{TRUE})?
#' @param map.unit \code{character}. Unit of position in \code{map}, i.e. 'cM'
#' for genetic distance or 'bp' for physical distance (default = 'cM').
#' @param repeated This column is used to identify the replications of the
#' phenotypic values. The unique values become the names of the third dimension
#' of the pheno object in the \code{gpData}. This argument is only required for
#' repeated measurements.
#' @param modCovar \code{vector} with \code{colnames} which identify columns
#' with covariables in \code{pheno}. This argument is only required for
#' repeated measurements.
#' @param na.string \code{character} or vector of \code{characters}. You can
#' specify values with which \code{NA} is coded in your geno object. In case
#' you read missing values from a file not as missing, but as character
#' strings. It can be specified more than one value for missings in a vector.
#' Default is \code{"NA"}.
#' @param cores \code{numeric}. Here you can specify the number of cores you
#' like to use.
#' @return Object of class \code{gpData} which is a \code{list} with the
#' following elements \item{covar}{\code{data.frame} with information on
#' individuals} \item{pheno}{\code{array} (individuals x traits x replications)
#' with phenotypic data} \item{geno}{\code{matrix} marker matrix containing
#' genotypic data. Columns (marker) are in the same order as in \code{map} (if
#' \code{reorderMap=TRUE}.) } \item{pedigree}{object of class \code{pedigree}}
#' \item{map}{\code{data.frame} with columns 'chr' and 'pos' and markers sorted
#' by 'pos' within 'chr'} \item{phenoCovars}{\code{array} with phenotypic
#' covariates} \item{info}{\code{list} with additional information on data
#' (coding of data, unit in \code{map}) From synbreed version 0.11-11 on the
#' function \code{\link{codeGeno}} adds here the package version which was used
#' to do the coding. There are differences in codings between version 0.10-11
#' and 0.11-0! }
#' @note In case of missing row names or column names in one item, information
#' is substituted from other elements (assuming the same order of
#' individuals/markers) and a warning specifying the assumptions is returned.
#' Please check them carefully.
#' @author Valentin Wimmer and Hans-Juergen Auinger with contributions be Peter
#' VandeHaar
#' @seealso \code{\link{codeGeno}}, \code{\link{summary.gpData}},
#' \code{\link{gpData2data.frame}}
#' @keywords manip
#' @examples
#'
#' set.seed(123)
#' # 9 plants with 2 traits
#' n <- 9 # only for n > 6
#' pheno <- data.frame(Yield = rnorm(n, 200, 5), Height = rnorm(n, 100, 1))
#' rownames(pheno) <- letters[1:n]
#'
#' # marker matrix
#' geno <- matrix(sample(c("AA", "AB", "BB", NA),
#'   size = n * 12, replace = TRUE,
#'   prob = c(0.6, 0.2, 0.1, 0.1)
#' ), nrow = n)
#' rownames(geno) <- letters[n:1]
#' colnames(geno) <- paste("M", 1:12, sep = "")
#'
#' # genetic map
#' # one SNP is not mapped (M5) and will therefore be removed
#' map <- data.frame(chr = rep(1:3, each = 4), pos = rep(1:12))
#' map <- map[-5, ]
#' rownames(map) <- paste("M", c(1:4, 6:12), sep = "")
#'
#' # simulate pedigree
#' ped <- simul.pedigree(3, c(3, 3, n - 6))
#'
#' # combine in one object
#' gp <- create.gpData(pheno, geno, map, ped)
#' summary(gp)
#'
#'
#' # 9 plants with 2 traits , 3 replications
#' n <- 9 #
#' pheno <- data.frame(
#'   ID = rep(letters[1:n], 3), rep = rep(1:3, each = n),
#'   Yield = rnorm(3 * n, 200, 5), Height = rnorm(3 * n, 100, 1)
#' )
#'
#' # combine in one object
#' gp2 <- create.gpData(pheno, geno, map, repeated = "rep")
#' summary(gp2)
#' @export create.gpData
#' @importFrom doBy orderBy
#' @importFrom methods is
#' @importFrom stats as.formula rnorm
#' @importFrom utils sessionInfo
#'
create.gpData <- function(pheno = NULL, geno = NULL, map = NULL, pedigree = NULL, family = NULL, covar = NULL,
                          reorderMap = TRUE, map.unit = "cM", repeated = NULL, modCovar = NULL, na.string = "NA", cores = 1) {
  infoCall <- match.call()
  # start with some checks on data
  # geno as matrix but not data.frame (storage)
  if (!map.unit %in% c("cM", "bp", "kb", "Mb")) warning("The measurement unit for the positions in the map should be either 'cM', 'bp', 'kb' or 'Mb'")
  if (!is.null(geno)) {
    if (is.data.frame(geno)) {
      geno <- as.matrix(geno)
      # if(any(duplicated(geno,MARGIN=1))) warning("individuals with duplicated genotypes")
    }
    geno[geno %in% na.string] <- NA
    if (!is.matrix(geno)) stop("geno must be a matrix or data.frame, not a ", class(geno))
    if (anyDuplicated(rownames(geno))) {
      stop(paste("In", substitute(geno), " are duplicated names of genotypes in rownames!"))
    }
    if (is.null(rownames(geno)) && is.null(pheno)) stop("rownames(geno) cannot be NULL unless a pheno of the same length is supplied")
    if (is.null(colnames(geno)) && is.null(map)) warning("colnames(geno) should not be NULL unless a map of the same length is supplied")
  }

  if (!is.null(map)) {
    if (!is.data.frame(map)) stop("map must be a data.frame, not a", class(map), "object")
    # as a data.frame, map already has rownames
    if (!all(c("chr", "pos") %in% colnames(map))) stop("colnames(map) must include 'chr' and 'pos'")
    # test if positions in map are numeric
    if (!is.numeric(map$pos)) stop("Position informations have to be numeric values!")

    # test if chr in map is numeric or character
    if (!(is.numeric(map$chr) | is.character(map$chr))) {
      warning(paste("Chromosome information (in map$chr) should be numeric or character, not", class(map$chr)))
    }
  } else {
    map.unit <- NULL
  }

  # match geno and map
  if (!is.null(geno) & !is.null(map)) {
    if (!all(colnames(geno) %in% rownames(map))) {
      warning("not all markers in 'geno' mapped in 'map'. gaps filled with 'NA' \n")
    }
    if (ncol(geno) == nrow(map) && is.null(colnames(geno))) {
      # assuming same markers in geno and map
      warning("missing colnames in 'geno': assuming to be identical as rownames in 'map' because of identical length \n")
      colnames(geno) <- rownames(map)
    } else {
      if (!identical(colnames(geno), rownames(map))) {
        # markers must be reordered
        map <- data.frame(
          chr = map$chr[match(colnames(geno), rownames(map))],
          pos = map$pos[match(colnames(geno), rownames(map))],
          row.names = colnames(geno)
        )
      }
    }
  }

  phenoCovars <- NULL
  attrModCovars <- NULL
  if (!is.null(pheno)) {
    if (2 != length(dim(pheno))) stop("pheno must be 2-dimensional, not ", length(dim(pheno)), "-dimensional")
    classList <- unlist(lapply(pheno, class))
    if (!all((classList[!names(classList) %in% repeated & !names(classList) %in% modCovar])[-1] %in% c("numeric", "integer"))) stop("Trait values have to be numeric!")
    # repeated measures? Use rownames of pheno as identifier for genotypes
    if (is.null(repeated)) {
      if (is.null(rownames(pheno))) stop("rownames(pheno) cannot be null when not using repeated measures!")
      if (anyDuplicated(rownames(pheno))) warning("rownames in pheno should be unique when not using repeated measures!")
      add <- 10^ceiling(log10(nrow(pheno)))
      # if the rownames are numbers, then temporarily add a big number to allow alphabetical sorting.
      if (all(rownames(pheno) %in% 1:nrow(pheno))) rownames(pheno) <- add + as.numeric(rownames(pheno)) else add <- NULL
      if (dim(pheno)[2] == 1) { # only a vector of traits
        phenoNames <- dimnames(pheno)
        arrPheno <- array(pheno[order(phenoNames[[1]]), ], dim = c(length(phenoNames[[1]]), 1, 1))
        dimnames(arrPheno) <- list(phenoNames[[1]][order(phenoNames[[1]])], phenoNames[[2]], "1")
      } else { # more than one trait, still unreplicated
        pheno <- pheno[order(rownames(pheno)), ]
        pheno.not.modCovar <- pheno[, !(colnames(pheno) %in% modCovar)]
        arrPheno <- array(as.matrix(pheno.not.modCovar), dim = c(dim(pheno.not.modCovar), 1))
        dimnames(arrPheno) <- list(rownames(pheno), colnames(pheno.not.modCovar), "1")
      }
      if (!is.null(add)) dimnames(arrPheno)[[1]] <- as.numeric(dimnames(arrPheno)[[1]]) - add
      if (!is.null(modCovar)) {
        arrModCovars <- array(1, dim = c(dim(pheno[, colnames(pheno) %in% modCovar]), 1))
        dimnames(arrModCovars) <- list(dimnames(arrPheno)[[1]], colnames(pheno)[colnames(pheno) %in% modCovar], "1")
        for (i in colnames(pheno)[colnames(pheno) %in% modCovar]) {
          arrModCovars[, i, 1] <- pheno[, i]
        }
      }
    } else { # a vector with replication identifier is applied. The first column is the identifier for genotypes
      dim3 <- data.frame(unique(pheno[, repeated]))
      colnames(dim3) <- repeated
      dim3 <- orderBy(as.formula(paste("~", paste(repeated, collapse = " + "))), data = dim3)
      for (i in 1:ncol(dim3)) dim3[, i] <- as.character(dim3[, i])
      rownam <- sort(unique(pheno[, 1]))
      if (!is.null(modCovar)) repeated <- unique(c(repeated, modCovar))
      arrPheno <- array(NA, dim = c(length(rownam), ncol(pheno) - (1 + length(repeated)), nrow(dim3)))
      dimnames(arrPheno) <- list(rownam, (colnames(pheno)[!colnames(pheno) %in% repeated])[-1], as.character(apply(dim3, 1, paste, collapse = "_")))
      for (i in 1:nrow(dim3)) {
        vec.bool <- apply(as.matrix(pheno[, colnames(dim3)]) == as.matrix(dim3[rep(i, nrow(pheno)), ]), 1, all)
        arrPheno[as.character(pheno[vec.bool, 1]), , i] <- as.matrix(pheno[vec.bool, (colnames(pheno)[!colnames(pheno) %in% repeated])[-1]])
      }
      if (!is.null(modCovar)) { # strip out covariates
        arrModCovars <- arrPheno[, rep(1, length(modCovar)), ]
        dimnames(arrModCovars)[[2]] <- colnames(pheno)[colnames(pheno) %in% modCovar]
        for (i in 1:nrow(dim3)) {
          vec.bool <- apply(matrix(pheno[, colnames(dim3)] == dim3[rep(i, nrow(pheno)), ], ncol = ncol(dim3)), 1, all)
          arrModCovars[as.character(pheno[vec.bool, 1]), , i] <- as.matrix(pheno[vec.bool, colnames(pheno)[colnames(pheno) %in% modCovar]])
        }
      }
    }
    if (!is.null(modCovar)) { # take the correct class of the covariates
      phenoCovars <- arrModCovars
      attrModCovars <- classList[dimnames(arrModCovars)[[2]]]
      for (i in names(attrModCovars)) {
        if (attrModCovars[i] != "numeric") {
          attrModCovars[i] <- "factor"
        }
      }
    }
    pheno <- arrPheno
    rm(arrPheno)
  }

  # match geno and pheno
  if (!is.null(geno) & !is.null(pheno)) {
    if (is.null(dimnames(pheno)[[1]]) | is.null(rownames(geno))) {
      if (dim(pheno)[1] == nrow(geno)) {
        warning("assuming identical order of genotypes in 'pheno' and 'geno' because of identical length.\nControll the Output! There is no warranty of correctness!\n")
        if (is.null(dimnames(pheno)[[1]])) {
          dimnames(pheno)[[1]] <- rownames(geno)
        } else {
          rownames(geno) <- dimnames(pheno)[[1]]
        }
        if (is.null(dimnames(pheno)[[1]]) & is.null(rownames(geno))) dimnames(pheno)[[1]] <- rownames(geno) <- paste0("ID", 10^ceiling(log10(nrow(geno))) + 1:nrow(geno))
      } else {
        stop("missing rownames (animal IDs) for either 'pheno' or 'geno' and lengths do not agree.")
      }
      # now geno and pheno have rownames
    }
  }
  # sort geno by rownames (alphabetical order)
  if (!is.null(geno)) {
    if (all(row.names(geno) %in% 1:nrow(geno))) {
      geno <- geno[order(as.numeric(row.names(geno))), ]
    } else {
      geno <- geno[order(row.names(geno)), ]
    }
  }

  # sort markers by chromosome and position within chromosome
  if (!is.null(map)) {
    if (reorderMap) {
      map$sor <- substr(map$chr, nchar(as.character(map$chr)), nchar(as.character(map$chr)))
      if (any(unique(map$sor)[!is.na(unique(map$sor))] %in% 0:9)) map$sor <- 1
      # first order by rownames in alphabetical order (important for SNPs with the same position)
      map <- map[order(as.character(rownames(map))), ]
      map <- orderBy(~ sor + chr + pos, data = map)
      map$sor <- NULL
      # sortcolumns in geno, too
      geno <- geno[, rownames(map)]
    }
    map$pos[is.na(map$chr)] <- NA
    class(map) <- c("GenMap", "data.frame")
  }

  if (!is.null(pedigree)) {
    if (!is(pedigree, "pedigree")) warning(paste("object pedigree should be of class 'pedigree', not class", paste(class(pedigree), collapse = " ")))
  }

  # return object
  obj <- list(covar = NULL, pheno = pheno, geno = geno, map = map, pedigree = pedigree, phenoCovars = phenoCovars)

  # add information to element covar
  # sort all available individuals
  ids <- sort(unique(c(dimnames(obj$pheno)[[1]], rownames(obj$geno), as.character(obj$pedigree$ID))))
  if (all(ids %in% 1:length(ids))) ids <- sort(as.numeric(ids))

  obj$covar <- data.frame(
    id = ids,
    phenotyped = ids %in% dimnames(obj$pheno)[[1]],
    genotyped = ids %in% rownames(obj$geno),
    stringsAsFactors = FALSE
  )

  # family information for genotyped indviduals
  if (!is.null(family)) {
    if (!is.data.frame(family) & !is.matrix(family)) stop("family must be either a data.frame or a matrix, not a ", class(family))
    colnames(family)[1] <- "family"
    family$id <- as.character(rownames(family))
    obj$covar <- merge(obj$covar, family, by = "id", all = TRUE)
    obj$covar$genotyped[is.na(obj$covar$genotyped)] <- FALSE
    obj$covar$phenotyped[is.na(obj$covar$phenotyped)] <- FALSE
  } else {
    obj$covar$family <- rep(NA, nrow(obj$covar))
  }

  # add covar from arguments, if available
  if (!is.null(covar)) {
    if (is.matrix(covar)) {
      if (is.null(rownames(covar))) warning("the supplied covar's rownames will default to 1:nrow(covar), which is likely not correct.  Inspect the resulting $covar.")
      covar <- as.data.frame(covar)
    }
    if (!is.data.frame(covar)) stop("covar must be a data.frame, not a ", class(covar))
    # do not use any existing columns named 'genotyped', 'phenotyped' or 'id'
    covar <- covar[!colnames(covar) %in% c("genotyped", "phenotyped", "id", "family")]
    # merge with existing data
    obj$covar <- merge(obj$covar, covar, by.x = 1, by.y = 0, all = TRUE)
  }

  if (!is.null(obj$pedigree)) {
    obj$covar <- obj$covar[match(obj$pedigree$ID, obj$covar$id), ]
  }

  # further information
  obj$info$map.unit <- map.unit
  obj$info$codeGeno <- FALSE
  obj$info$attrPhenoCovars <- attrModCovars
  obj$info$version <- paste("gpData object was created by synbreed version", sessionInfo()$otherPkgs$synbreed$Version)
  obj$info$Call <- infoCall

  # return object of class 'gpData'
  class(obj) <- "gpData"
  return(obj)
}

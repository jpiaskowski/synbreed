
#' synbreed package: A Framework for the Analysis of Genomic
#' Prediction Data Using R
#'
#' The data flow in synbreed is guided by a single, unified data object of class
#' \code{gpData} (‘genomic prediction Data’) which is used for storage of multiple data
#' sources. This includes an array for the phenotypes (individual × trait ×
#' replication) and a matrix for the marker genotypes (individual × marker
#' scores). If required, this structure can be extended to include pedigree
#' information and a marker map. All analysis functions are based on this data
#' structure. A key feature of the synbreed package is the generality of the
#' class \code{gpData} which is suitable for a wide range of statistical methods using
#' genotypic and phenotypic data. Moreover, it is very convenient to store and
#' share objects of class gpData. Any object of class gpData can be converted to
#' class \code{cross} in the R package 'qtl' (Broman et al., 2003) and vice versa or to a
#' \code{data.frame).
#'
#' Main Functionality:
#'
#' - Data Preparation
#' - Visualisation
#' - `create.gpData`
#' - `gp.Mod`
#' - Genetic Modelling
#' - and more
#'
#' @name synbreed-package
#' @docType package
#' @keywords package
#'
#' @examples
#' # (better put something here)
NULL


#' Dairy cattle data
#'
#' Data set contains genotypic, phenotypic, map and pedigree data of 500 bulls.
#' All individuals are labeled with an unique ID, starting with ID1430 and
#' ending with ID1929. Genotypic and pedigree data are from a real cattle
#' data set while phenotypes were simulated.  Pedigree information is
#' available for parents and grandparents of the phenotyped
#' individuals.
#'
#' There are two quantitative phenotypes in this data set. The heritabilities of these
#' traits, 0.41 and 0.66, were estimated with a pedigree-based animal model using
#' the data on hand.
#'
#' Genotypic data consists of 7250 biallelic SNP markers for every phenotyped
#' individual with missing data included. SNPs are mapped across all 29
#' autosomes. Distances in the SNP map are given in megabases (Mb).
#'
#'
#' @name cattle
#' @docType data
#' @format object of class \code{gpData}
#' @keywords datasets
#' @examples
#'
#' \dontrun{
#' library(synbreed)
#' data(cattle)
#' summary(cattle)
#' }
#'
NULL



#' Simulated maize data
#'
#' This is a simulated dataset of a maize breeding program. Data comprise 1250
#' doubled haploid lines that were genotyped with 1117 polymorphic SNP
#' markers and phenotyped in a testcross with a single tester for one
#' quantitative trait. All individuals are labeled with a unique ID, starting
#' from 11360 to 12609. Markers are distributed along all 10 chromosomes of
#' maize. Pedigree information starts with base population and is available up
#' to 15 generations. The 1250 lines belong to 25 full sib families with 50
#' individuals in each family. In the simulation of true breeding values (TBV),
#' 1000 biallelic quantitative trait loci (QTL) with equal and additive (no
#' dominance or epistasis) effects were generated. True breeding values for
#' individuals were calculated according to \deqn{tbv=\sum_{k=1}^{1000}
#' QTL_k}{TBV=\sum QTL(k)} where \eqn{QTL_k}{QTL(k)} is the effect of the
#' \eqn{k}-th QTL. Phenotypic values were simulated according to
#' \deqn{y_i=tbv_i + \epsilon_i}{trait=tbv+e} where \eqn{\epsilon_i \sim
#' N(0,\sigma^2)}{e = N(0,sigma2)}. The value for \eqn{\sigma^2}{sigma2} was
#' chosen in a way that a given plot heritability of \eqn{h^2=0.197}{h2=0.197}
#' is realized. Note that true breeding values for 1250 phenotyped lines are
#' stored as \code{tbv} in \code{covar} of \code{gpData} object. Reported
#' phenotypic values of lines are adjusted values testcross means for yield
#' [dt/ha] evaluated in 3 locations.
#'
#'
#' @name maize
#' @docType data
#' @format object of class \code{gpData}
#' @keywords datasets
#' @examples
#'
#' \dontrun{
#' library(synbreed)
#' data(maize)
#' summary(maize)
#' }
#'
NULL



#' Heterogenous stock mice population
#'
#' Data set comprises public available data of 2527 (1293 males and 1234
#' females) heterogenous stock mice derived from eight inbred strains (A/J,
#' AKR/J, BALBc/J, CBA/J, C3H/HeJ, C57BL/6J, DBA/2J and LP/J) followed by 50
#' generations of pseudorandom mating. All individuals are labeled with a
#' unique ID, starting with \code{A048005080}. For all individuals, family, sex
#' (females=0, males=1), month of birth (1-12), birth year, coat color, cage
#' density and litter is available and stored in \code{covar}.
#'
#' The measured traits are described in Solberg et al. (2006). Here, the body
#' weight at age of 6 weeks [g] and growth slope between 6 and 10 weeks age
#' [g/day] are available. The heritabilities of these traits are reported as
#' 0.74 and 0.30, respectively (Valdar et al, 2006b). Phenotypic data was taken
#' from \url{http://mus.well.ox.ac.uk/GSCAN/HS_PHENOTYPES/Weight.txt}.
#'
#' Genotypic data consists of 12545 biallelic SNP markers and is available for
#' 1940 individuals. Raw genotypic data from
#' \url{http://mus.well.ox.ac.uk/GSCAN/HS_GENOTYPES/} is given in the
#' \code{Ped-File Format} with two columns for each marker. Both alleles were
#' combined to a single genotype for each marker in \code{mice} data. The SNPs
#' are mapped in a sex-averaged genetic map with distances given in centimorgan
#' (Shifman et al. (2006)). SNPs are mapped across all 19 autosomes and
#' X-chromosome where distances between adjacent markers vary form 0 to 3 cM.
#'
#'
#' @name mice
#' @docType data
#' @format object of class \code{gpData}
#' @references
#'
#' Shifman S, Bell JT, Copley RR, Taylor MS, Williams RW, et al.
#' (2006) A High-Resolution Single Nucleotide Polymorphism Genetic Map of the
#' Mouse Genome. PLoS Biol 4(12)
#'
#' Solberg L.C. et al. (2006), A protocol for high-throughput phenotyping,
#' suitable for quantitative trait analysis in mice. Mamm. Genome 17, 129-146
#'
#' Valdar W, Solberg LC, Gauguier D, Burnett S, Klenerman P, Cookson WO, Taylor
#' MS, Rawlins JN, Mott R, Flint J. (2006a) Genome-wide genetic association of
#' complex traits in heterogeneous stock mice. Nat Genet. 8, 879-887.
#'
#' Valdar W, Solberg LC, Gauguier D, Cookson WO, Rawlins NJ, Mott R, Flint
#' J.(2006b) Genetic and environmental effects on complex traits in mice.
#' Genetics 175, 959-984
#'
#' @source Welcome Trust Centre for Human Genetics, Oxford University, data
#' available from \url{http://gscan.well.ox.ac.uk}
#' @keywords datasets
#' @examples
#'
#' \dontrun{
#' library(synbreed)
#' data(mice)
#' summary(mice)
#' }
#'
NULL

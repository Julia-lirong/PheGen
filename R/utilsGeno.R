


#' Compute allele frequencies from genotype data.
#'
#' @rdname utilizeGeno
#' @param geno [N x M] Vector of length [N] samples with genotypes of [M] single
#'   bi-allelic genetic variant/SNP encoded as 0,1 and 2.
#' @return Vector with minor allele frequencies.
#' @importFrom Hmisc impute
#' @examples
#' geno <- simulateGeno(N=10, NrSNP=100, frequency= c(0.3,0.5))$genotypes
#' allelefreq <- getmaf(geno)

getmaf <- function(geno) {
  if (any(is.na(geno))) {
    stop("Missing genotypes found, remove missing genotypes or impute by mean")
  }
  if (any(is.na(geno))) {
    geno <- round(apply(as.matrix(geno), 2, impute, fun = mean))
  }
  p <- sapply(1:ncol(geno), function(x) {
    sum(geno[, x])/(2 * nrow(geno))
  })
  p <- sapply(1:length(p), function(x) {
    min(p[x], 1 - p[x])
  })
  return(p)
}


#' Standardise genotypes.
#'
#' Genotypes are standardised as described in Yang et al:
#' snp_standardised = (snp - 2 * ref_allele_freq)/
#' sqrt(2 * ref_allele_freq * alt_allele_freq).
#'
#' Missing genotypes can be mean-imputed and rounded to nearest integer
#' before standardisation. If genotypes contain missing values and impute is set
#' to FALSE, \code{standardiseGenotypes} will return an error.
#'
#' @rdname utilizeGeno
#' @param geno [N x NrSNP] Matrix/dataframe of genotypes [integer]/[double].
#' @param impute [logical] Indicating if missing genotypes should be imputed; if
#'   set FALSE and data contains missing values,  \code{standardiseGenotypes} will
#'   return an error.
#' @return [N x NrSNP] Matrix of standardised genotypes [double].
#' @references Yang, J., Lee, S.H., Goddard, M.E., Visscher, P.M. (2011) GCTA:
#'   a tool for genome-wide complex trait analysis, AJHG: 88
#' @importFrom Hmisc impute
#' @examples
#' geno <- simulateGeno(N=10, NrSNP=100, frequency= c(0.3,0.5))$genotypes
#' geno_sd <- standardiseGeno(geno)
standardiseGeno <- function(geno, impute = TRUE) {
  if (any(is.na(geno)) & !impute) {
    stop("Missing genotypes found and impute=FALSE, cannot standardise",
         "genotypes; remove missing genotypes or set impute=TRUE for mean",
         "imputation of genotypes")
  }
  if (any(is.na(geno)) & impute) {
    geno <- round(apply(as.matrix(geno), 2, impute, fun = mean))
  }
  p <- getmaf(geno)
  geno <- sapply(1:ncol(geno), function(i) {
    if (p[i] != 0) {
      (geno[, i] - 2 * p[i])/sqrt(2 * p[i] * (1 - p[i]))
    } else {
      geno[, i]
    }
  })

  return(geno)
}








############################################################
######### function from SKAT pacakge #######################
############################################################

# delete the snp with maf=0
Delet_maf_hap <- function(Haplotype, SNPInfo) {
  Marker.MAF.ALL <- colMeans(Haplotype)
  # delete the variant with maf=0
  id.non <- which(Marker.MAF.ALL == 0)
  if (length(id.non) > 0) {
    Haplotype <- Haplotype[, -id.non]
    SNPInfo <- SNPInfo[-id.non, ]
  }
  return(list(Haplotype = Haplotype, SNPInfo = SNPInfo))
}


# SNP.Dist: SNP.Location LIST,
# SubRegion.Length: SubRegion.Length, correspongding to the position, then choose the varaince in the base pair
Get_RandomRegion <- function(SNP.Dist, SubRegion.Length) {
  if (SubRegion.Length < 0) {
    return(1:length(SNP.Dist))
  }
  total.first <- min(SNP.Dist)
  total.last <- max(SNP.Dist)
  total.length <- total.last - total.first
  Region.Start <- runif(1) * (total.length - SubRegion.Length) + total.first
  Region.End <- Region.Start + SubRegion.Length
  Marker.Idx1 <- which(SNP.Dist >= Region.Start)
  Marker.Idx2 <- which(SNP.Dist <= Region.End)
  IDX <- sort(intersect(Marker.Idx1, Marker.Idx2))
  return(IDX)
}






Get_Gene <- function(Haplotype, N, IDX.Marker) {
  n1 <- dim(Haplotype)[1]  # the number of Haplotype
  # p <- colMeans(Haplotype)
  H1 <- sample(1:n1, N, replace = TRUE)
  H2 <- sample(1:n1, N, replace = TRUE)
  gene <- Haplotype[H1, IDX.Marker] + Haplotype[H2, IDX.Marker]
  if (N==1) gene <- as.data.frame(t(gene))
  # set at least one genotype is nonzero
  ind <- which(colMeans(gene)/2 == 0)
  rand <- sample(1:N, length(ind), replace = TRUE)
  if (length(ind) > 0) {
    for (i in 1:length(ind)) {
      gene[rand[i], ind[i]] <- 1
    }
  }
  return(gene)
}




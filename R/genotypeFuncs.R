

#' Simulate genotypes based on predefined frequencies
#'
#' \code{simulateGeno()} Simulate genotypes in a region based on given minor allele frequency
#'
#' @aliases simulateGeno simulateHap simulateGeneHap
#' @rdname GenotypeSimulate
#' @param N Number of samples for which to simulate bi-allelic genotypes.
#' @param NrSNP Number of SNPs to simulate.
#' @param SNPfrequency Vector of allele frequencies [double] from which to
#'   sample.
#' @param sampleID Prefix [string] for naming samples (will be followed by
#'   sample number from 1 to N when constructing id_samples).
#' @param snpID Prefix [string] for naming SNPs (will be followed by SNP number
#'   from 1 to NrSNP when constructing id_snps).
#' @param is.standardise Logical value. If TRUE, standardize the genotype matrix.
#'   Default is FALSE.
#' @param verbose [boolean] If TRUE, progress info is printed to standard out.
#' @param \dots Arguments passed to the internal function.
#' @return Named list with [N x NrSNP] matrix of simulated genotypes
#'   (genotypes), their SNP frequencies (freq), a vector of sample IDs
#'   (id_samples) and a vector of SNP IDs (id_snps).
#' @import stats
#' @importFrom utils data
#' @export
#' @examples
#' data(frequencies)
#' snp <- simulateGeno(N=10, NrSNP=10)

simulateGeno <- function(N, NrSNP = NULL,
                         SNPfrequency = NULL, sampleID = "ID_",
                         snpID = "SNP_", is.standardise = FALSE,
                         verbose = TRUE, ...) {
  if (is.null(SNPfrequency)) {
    # vmessage(paste0("Using the default frequencies data"))
    data(frequencies, package = "PheGen")
  }
  SNPfrequency <- frequencies$MAF
  samples <- paste(sampleID, 1:N, sep = "")
  snps <- paste(snpID, 1:NrSNP, sep = "")
  # vmessage(c("Simulate", NrSNP, "SNPs for", N, "individuals"), verbose = verbose)
  freq <- sample(SNPfrequency, NrSNP, replace = TRUE)
  if (N > 0) {
    gene <- sapply(1:NrSNP, function(x) rbinom(N, 2, freq[x]))
    if (N == 1) gene <- as.data.frame(t(gene))
  } else {
    gene <- matrix(0, nrow = 0, ncol = NrSNP)
  }
  if (is.standardise) {
    gene <- standardiseGeno(geno = gene, ...)
  }
  SNPInfo <- data.frame(SNP = c(1:NrSNP), FREQ1 = freq)
  return(list(genotypes = gene, SNPInfo = SNPInfo))
}






#' Simulate genotypes in a region from the haplotype pool
#' @rdname GenotypeSimulate
#' @param Haplotype Haplotype pool (two objects, first is the hoplotypes,
#'   second object is the SNP information)
#' @importFrom utils data
#' @export
#'
#' @examples
#' gene <- simulateHap(N=100, NrSNP=20)

simulateHap <- function(N, NrSNP,
                        Haplotype = NULL,
                        is.standardise = FALSE, ...) {
  if (is.null(Haplotype)) {
    # vmessage(paste0("Using the default hoplotype data"))
    data(haplotypes, package = "PheGen")
  }
  Hap <- haplotypes$Haplotype
  IDX.Marker <- sample(1:ncol(Hap), NrSNP, replace = TRUE)
  H1 <- sample(1:nrow(Hap), N, replace = TRUE)
  H2 <- sample(1:nrow(Hap), N, replace = TRUE)
  gene <- Hap[H1, IDX.Marker] + Hap[H2, IDX.Marker]
  freq = colMeans(gene)/2
  if (is.standardise) {
    gene <- standardiseGeno(geno = gene, ...)
  }
  SNPInfo <- data.frame(SNP = c(1:NrSNP), FREQ1 = freq)

  return(list(genotypes = gene, SNPInfo = SNPInfo))

}




#' Simulate genotypes in a region from the haplotype pool
#' @rdname GenotypeSimulate
#' @param Haplotype Haplotype pool (two objects, first is the hoplotypes,
#'   second object is the SNP information)
#' @param SubRegion.Length A value of the length of subregions.
#'   If SubRegion.Length=-1 (default), the length of the subregion
#'   will be the same as the length of the whole region, so there will
#'   no random selection of subregions. This parameter is used in the gene-based
#'   assocaition.
#' @importFrom utils data
#' @export
#'
#' @examples
#' gene <- simulateGeneHap(N=100)

simulateGeneHap <- function(N, Haplotype = NULL,
                            SubRegion.Length = -1,
                            is.standardise = FALSE, ...) {
  if (is.null(Haplotype)) {
    # vmessage(paste0("Using the default hoplotype data"))
    data(haplotypes, package = "PheGen")
  }
  Hap <- haplotypes
  Haplotype <- Hap$Haplotype
  SNPInfo <- Hap$SNPInfo
  if (is.null(SubRegion.Length)) {
    stop("SubRegion.Length is a value between 79 to 199,956 ")
  } else if (SubRegion.Length == -1) {
    IDX.Marker <- c(1:ncol(Haplotype))
  } else {
    IDX.Marker <- Get_RandomRegion(SNP.Dist = SNPInfo$CHROM_POS, SubRegion.Length)
  }
  gene_infor <- SNPInfo[IDX.Marker, ]
  M <- length(IDX.Marker)
  gene <- Get_Gene(Haplotype = Haplotype, N, IDX.Marker)
  if (is.standardise) {
    gene <- standardiseGeno(geno = gene, ...)
  }
  return(list(genotypes = gene, SNPInfo = gene_infor))

}




#' Phenotype simulation.
#'
#' Phenotypes are modeled under a linear additive model where Y = GB + E,
#' GB the genetic variant effects, E is noise term.
#'
#' @rdname PhenotypeSimulate
#' @aliases PheSimulator
#' @param K Number [integer] of phenotypes to simulate.
#' @param N Vector of sample sizes for K phenotypes. If input is a single integer
#'   which means all phenotypes have same sample size.
#' @param NrSNP Number [integer] of genetic variants to simulate.
#' @param Ns The number of overlapped samples among these phenotypes.
#' @param genoMethod Name of methods to generate Genotypes, which are "SNPfrequency",
#'   , "Haplotype" or "Region". Default is "SNPfrequency".
#' @param SNPfrequency Vector of allele frequencies [double] from which to sample.
#' @param cSNP Proportion [double] of causal SNPs; used as genetic variant effects.
#'   default is 0.01.
#' @param pcorr [K x K] matrix of correlation [double] between phenotypes. If input
#'   is a single value which means the correlation between any pair of
#'   phenotypes are same. Default is 0.5.
#' @param pgenVar [K x K] matrix. Proportion of phenotype correlation influenced 
#'   by genetic variants. Default is 0.6. It is used with \code{pcorr} parameters.
#' @param h2s K by 1 vector of genomic heritability of different phenotypes.
#'   Each element means the proportion of variance of a phenotype that can be
#'   explained by a linear regression on a set of variants. Default is 0.5.
#' @param theta Proportion [double] of variance of shared genetic variant effects.
#' @param weight [NrSNP x NrSNP] diagonal matrix of the weight which is 
#'   assigned to the variants. Default is 1.
#' @param Haplotype Haplotype pool (two objects, first is the hoplotypes,
#'   second object is the SNP information)
#' @param method String specifying the matrix decomposition used to determine
#'   the matrix root of sigma. Possible methods are eigenvalue decomposition
#'   ("eigen", default), singular value decomposition ("svd"), and
#'   Cholesky decomposition ("chol"). The Cholesky is typically fastest,
#'   not by much though.
#' @param \dots Arguments passed to the internal function.
#'
#' @return Named list of i) the final simulated phenotype components
#'   (phenoComponentsFinal), ii) a named list of raw components (rawComponents)
#'   used for genetic effect simulation (genotypes and/or kinship, eigenvalues
#'   and eigenvectors of kinship)
#' @export
#'
#' @examples
#' K <- 5; N<- c(100,200,300,150,200); Ns<-100; NrSNP=1000
#' h2s <- rep(0.6,K); pcorr <- 0.5
#' out <- PheSimulator(K, N, NrSNP, Ns, pcorr=pcorr, h2s=h2s, genoMethod = "Haplotype")


PheSimulator <- function(K, N, NrSNP, Ns, genoMethod = "SNPfrequency",
                         SNPfrequency = NULL, cSNP = 0.02, 
                         pcorr = NULL, pgenVar = 0.6, h2s = NULL,
                         theta = 0.8, weight = 1, Haplotype, method = "eigen", ...) {
  
  if (is.null(pcorr) | is.null(h2s)) {
    stop(paste("Phenotypic correlation or heritability is missing"))
  }
  
  # Generate genotype
  if (genoMethod == "SNPfrequency") {
    if (Ns > 0) {
      gene <- simulateGeno(N = Ns, NrSNP, is.standardise = FALSE, ...)$genotypes
      genes <- sapply(1:K, function(i) {
        rbind(gene, simulateGeno(N = N[i] - Ns, NrSNP, is.standardise = FALSE, ...)$genotypes)
      })
    }
    Haps <- NULL
  } else if (genoMethod == "Region") {
    gene <- simulateGeneHap(N = Ns, SubRegion.Length = 10, is.standardise = FALSE, ...)$genotypes
    genes <- sapply(1:K, function(i) {
      rbind(gene, simulateGeneHap(N = N[i] - Ns, SubRegion.Length = 10, is.standardise = FALSE, ...)$genotypes)
    })
    Haps <- NULL
  } else if (genoMethod == "Haplotype") {
    temp <- simulateHap(N = Ns, NrSNP, is.standardise = FALSE, ...)
    gene <- temp$genotypes
    Haptemp <- temp$Haplotype
    genes <- Haps <- vector(mode='list', length=K)
    for (i in 1:K) {
      temp <- simulateHap(N = N[i] - Ns, NrSNP, is.standardise = FALSE, ...)
      genes[[i]] <- rbind(gene, temp$genotypes)
      Haps[[i]] <- rbind(Haptemp, temp$Haplotype)
    }
  }
  
  gene_st <- lapply(1:K, function(i){
    standardiseGeno(geno = genes[[i]])
  })
  
  # generate effect size
  M <- ceiling(NrSNP * cSNP)
  id.causal <- sample(1:NrSNP, M, replace = FALSE)
  if (is.null(weight)) {
    W <- diag(rep(1, M))
  } else if (length(weight) == 1) {
    W <- diag(weight, nrow = M, ncol = M)
  }
  
  if (length(pcorr) == 1) {
    V.a <- matrix(0, nrow = K, ncol = K)
    upper_tri <- upper.tri(V.a, diag = FALSE)
    V.a[upper_tri] <- cumprod(rep(pcorr, K * (K - 1)/2))
    # V.a[upper_tri] <- rep(pcorr, K* (K - 1)/2)
    V.a <- V.a + t(V.a)
    diag(V.a) <- 1
  }
  if (length(pgenVar) == 1) {
    pgenVar <- matrix(pgenVar, nrow = K, ncol = K)
  }
  V.g <- V.a * pgenVar/M
  diag(V.g) <- h2s/M
  V.e <- V.a - V.g*M
  
  B <- matrix(0, nrow = NrSNP, ncol = K)
  B[id.causal, ] <- geneticEffect(M, K, W, V.g, method)
  
  Nmax <- max(N)
  E <- noiseEffect(Nmax, K, V.e, method)
  
  Y <- lapply(1:K, function(i) {
    gene_st[[i]] %*% B[, i] + E[1:N[i], i]
  })
  Y <- lapply(1:K, function(i){
    scale(Y[[i]])
  })
  return(list(Genotype = genes, Phenotype = Y, Hap = Haps))
}



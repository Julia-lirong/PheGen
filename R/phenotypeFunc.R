
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
#'   or "Haplotype". Default is "SNPfrequency".
#' @param SNPfrequency Vector of allele frequencies [double] from which to sample.
#' @param cSNP Proportion [double] of causal SNPs; used as genetic variant effects.
#'   default is 0.01.
#' @param genVar [K x K] matrix of genetic covariance between any pair of phenotypes.
#' @param noiseVar [K x K] matrix of noise covaraince between any pair od phenotypes.
#' @param pcorr [K x K] matrix of correlation [double] between phenotypes. If input
#'   is a single value which means the correlation between any pair of
#'   phenotypes are same. Default is 0.5.
#' @param pgenVar Proportion of phenotype correlation influenced by genetic
#'   variants. Default is 0.6. It is used with \code{pcorr} parameters.
#' @param h2s K by 1 vector of genomic heritability of different phenotypes.
#'   Each element means the proportion of variance of a phenotype that can be
#'   explained by a linear regression on a set of variants. Default is 0.5.
#' @param theta Proportion [double] of variance of shared genetic variant effects.
#' @param pTraits proportion [double] of traits which are correlated to each other.
#'   Allows to simulate for instance different levels of pleiotropy. Using
#'   ceiling(K*pTraits) as the number of phenotypes correlated to each other.
#' @param Haplotype Haplotype pool (two objects, first is the hoplotypes,
#'   second object is the SNP information)
#'
#' @return Named list of i) the final simulated phenotype components
#'   (phenoComponentsFinal), ii) a named list of raw components (rawComponents)
#'   used for genetic effect simulation (genotypes and/or kinship, eigenvalues
#'   and eigenvectors of kinship)
#' @export
#'
#' @examples
#' K <- 5; N<- c(100,200,300,150,200); M<-160; Ns<-100
#'
#'
PheSimulator <- function(K,N,NrSNP,Ns, genoMethod="SNPfrequency",
                         SNPfrequency = NULL,cSNP=0.01,
                         genVar=NULL, noiseVar=Null,
                         pcorr=NULL, pgenVar= 0.6,h2s=NULL,
                         theta=0.8,
                         pTraits=1,Haplotype,...){

  if (is.null(genVar)) {
    if (is.null(noiseVar) & is.null(pcorr)) {
      stop(paste("Neither genVar nor noiseVar are provided, thus",
                 "proportion of variance from genetics cannot be deduced"))
    }
  }
  if (is.null(noiseVar)) {
    noiseVar <- pcorr - genVar
    if (min(noiseVar)<0){
      stop(paste("noiseVar can not greater than 1"))
    }
  }

  if (noiseVar + genVar > 1) {
    stop("Sum of genetic variance and noise variance is greater than 1")
  }

  # Generate genotype
  if (genoMethod=="SNPfrequency") {
    if (Ns > 0) {
      gene <- simulateGeno(N=Ns, NrSNP)$genotypes
      genes <- sapply(1:K, function(i){
        rbind(gene,simulateGeno(N=N[i]-Ns, NrSNP)$genotypes)
      })
    }
  } else if (genoMethod=="Haplotype"){
    gene <- simulateGeneHap(N=Ns,SubRegion.Length = 10)$genotypes
    genes <- sapply(1:K, function(i){
      rbind(gene,simulateGeneHap(N=N[i]-Ns, SubRegion.Length = 10)$genotypes)
    })
  }




  corr_vec <- cumprod(rep(pcorr, P - 1))
  N <- 10
  M <- 20
  K <- 5
  M.causal <- 15
  W <- diag(rep(1,M))
  mu <- rep(0, M)
  V.g <- V.e <-  matrix(0,nrow = K,ncol = K)
  #h2 <- runif(K,0,1)
  h2 <- rep(0.9, K)
  #phenotype correlation matrix due to genetic effect
  upper_tri <- upper.tri(V.g, diag=FALSE)
  V.g[upper_tri] <- runif(sum(rep(1:(K-1))),0,1)
  # V.g[upper_tri] <- 0.6
  V.g <- V.g + t(V.g)
  diag(V.g) <- 1
  V.g <- diag(sqrt(h2)) %*% V.g %*% diag(sqrt(h2))
  V.g <- V.g/M.causal

  ####
  upper_tri <- upper.tri(V.e, diag=FALSE)
  V.e[upper_tri] <- 0.9
  V.e <- V.e + t(V.e)
  diag(V.e) <- 1
  V.e <- diag(sqrt(1-h2)) %*% V.e %*% diag(sqrt(1-h2))

  V.g + V.e
}



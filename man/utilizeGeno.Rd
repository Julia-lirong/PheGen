% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilsGeno.R
\name{getmaf}
\alias{getmaf}
\alias{standardiseGeno}
\title{Compute allele frequencies from genotype data.}
\usage{
getmaf(geno)

standardiseGeno(geno, impute = TRUE)
}
\arguments{
\item{geno}{[N x NrSNP] Matrix/dataframe of genotypes [integer]/[double].}

\item{impute}{[logical] Indicating if missing genotypes should be imputed; if
set FALSE and data contains missing values,  \code{standardiseGenotypes} will
return an error.}
}
\value{
Vector with minor allele frequencies.

[N x NrSNP] Matrix of standardised genotypes [double].
}
\description{
Genotypes are standardised as described in Yang et al:
snp_standardised = (snp - 2 * ref_allele_freq)/
sqrt(2 * ref_allele_freq * alt_allele_freq).
}
\details{
Missing genotypes can be mean-imputed and rounded to nearest integer
before standardisation. If genotypes contain missing values and impute is set
to FALSE, \code{standardiseGenotypes} will return an error.
}
\examples{
#geno <- simulateGeno(N=10, NrSNP=100, frequency= c(0.3,0.5))$genotypes
#allelefreq <- getmaf(geno)
#geno <- simulateGeno(N=10, NrSNP=100, frequency= c(0.3,0.5))$genotypes
#geno_sd <- standardiseGeno(geno)
}
\references{
Yang, J., Lee, S.H., Goddard, M.E., Visscher, P.M. (2011) GCTA:
  a tool for genome-wide complex trait analysis, AJHG: 88
}

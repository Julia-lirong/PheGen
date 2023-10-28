


#' GWAS summary statistics calculation
#'
#' @param X Genotypes dataframe list with K objects
#' @param Y Phenotype dataframe list with K objects
#' @return GWAS summary statistics dataframe
#' @export
#'
#' @examples
#' K <- 5; N<- c(100,200,300,150,200); Ns<-100; NrSNP=1000
#' h2s <- rep(0.6,K); pcorr <- 0.5
#' out <- PheSimulator(K, N, NrSNP, Ns, pcorr=pcorr, h2s=h2s)
#' X <- out$Genotype
#' Y <- out$Phenotype
#' out <- GWASfunc(X,Y)


GWASfunc <- function(X, Y) {
  K <- length(Y)
  res <- lapply(1:K, function(i.K) {
    SCORE(x = X[[i.K]], y = Y[[i.K]])
  })
  return(res)
}







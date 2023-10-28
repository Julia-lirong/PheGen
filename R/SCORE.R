

#' SCORE function
#'
#' Score test for testing the association between phenotypes and SNPs.
#' Note that the phenotypes could be the binary phenotypes or
#' quantitative phenotypes.
#'
#' \code{SCORE} Default method;
#'
#' \code{SPA_SCORE} Score test by saddlepoint approximation for testing the
#'  association between only binary phenotypes and SNPs. Note that SPA_SCORE
#'  function can use to binary phenotypes with the unbalanced or extremely
#'  unbalanced cast-control ratios.
#'
#' @aliases SCORE SPA_SCORE Get_GWAS
#' @rdname SCORE
#' @param x The individual-level genotyes with dimension n individuals by M SNPs.
#'  It could be only one SNPs with a n by 1 genotype vector.
#' @param y The individual-level phentypes with dimesion n individuals by 1
#'  phenotypes.
#' @param cov The covariate matrix with dimesion n individuals by C covariates.
#'  defalut: there is no covariates.
#' @import stats
#' @return a dataframe of GWAS summary
#' @examples
#' #N <- 100; M <- 200; K <- 1
#' #x <- replicate(M, rbinom(N,2,0.3))
#' #y1 <- replicate(1, sample(c(0,1), N, replace = TRUE, prob = c(0.9, 0.1)))
#' #y2 <- replicate(1, rnorm(N))
#' #res <- SCORE(x, y2)
#' #res1 <- SPA_SCORE(x, y1)
#' #res2 <- Get_GWAS(x,y2)


SCORE <- function(x, y, cov = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  # Check dimension of the individual-level data: x and y
  if (nrow(x) != nrow(y)) {
    stop("Error: Please check the sample size of x and y.
         They should be the same!")
  }
  
  # Check the existing of the covariance matrix
  if (!is.null(cov)) {
    cov <- as.matrix(cov)
    if (nrow(cov) != nrow(x)) {
      stop("Error: Please check the sample size of covariance matrix,
           which should be the same as the sample size of  individual-level
           genotyps and phenotypes!")
    }
    x1 <- apply(x, 2, function(t) return(residuals(lm(t ~ 1 + cov))))
    y1 <- apply(y, 2, function(t) {
      if (all(t %in% c(0, 1))) {
        return(residuals(glm(t ~ 1 + cov, family = binomial(link = "logit"))))
      } else {
        return(residuals(lm(t ~ 1 + cov)))
      }
    })
    x <- as.matrix(x1)
    y <- as.matrix(y1)
  }
  
  n <- nrow(x)
  Tstat <- sqrt(n) * cor(y, x)
  pv <- pchisq(Tstat^2, 1, lower.tail = FALSE)
  out <- data.frame(pvalue = pv[1, ], Z = Tstat[1, ])
  return(out)
}



#' For binary phenotype
#' @rdname SCORE
#' @importFrom SPAtest ScoreTest_SPA
SPA_SCORE <- function(x, y, cov = NULL) {
  requireNamespace("SPAtest", quietly = TRUE)
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  # Check dimension of the individual-level data: x and y
  if (nrow(x) != nrow(y)) {
    stop("Error: Please check the sample size of x and y. They should be the same!")
  }
  # Check if binary phenotype
  if (!all(y %in% c(0, 1))) {
    stop("Error: SPA_SCORE only can be applied to the binary phenotypes")
  }
  # Check the existing of the covariance matrix
  if (!is.null(cov)) {
    cov <- as.matrix(cov)
    if (nrow(cov) != nrow(x)) {
      stop("Error: Please check the sample size of covariance matrix,
           which should be the same as the sample size of  individual-level
           genotyps and phenotypes!")
    }
  }
  M <- ncol(x)
  beta <- sd <- pvalue <- Z <- rep(NA, M)
  for (i.M in 1:M) {
    res <- ScoreTest_SPA(genos = x[, i.M], pheno = as.numeric(y[, 1]),
                         cov = cov, method = "fastSPA", beta.out = TRUE, beta.Cutoff = 1)
    pvalue[i.M] <- res$p.value
    beta[i.M] <- res$beta
    sd[i.M] <- res$SEbeta
    Z[i.M] <- res$beta/res$SEbeta
  }
  out <- data.frame(beta = beta, sd = sd, pvalue = pvalue, Z = Z)
  return(out)
}



#' For continuous phenotypes
#' @rdname SCORE
## get GWAS summary statistic
Get_GWAS <- function(x, y) {
  
  M <- ncol(x)
  beta <- sd <- pvalue <- Z <- rep(NA, M)
  for (i.M in 1:M) {
    model <- lm(y ~ x[, i.M])
    b <- unname(model$coefficients[2])  #effect size
    # Residual Standard error (Like Standard Deviation)
    k <- length(model$coefficients) - 1  #Subtract one to ignore intercept
    SSE <- sum(model$residuals^2)
    n <- length(model$residuals)
    SS_x <- (n - 1) * var(x[, i.M])
    se <- sqrt(SSE/((n - (1 + k)) * SS_x))  #Residual Standard Error
    t <- b/se
    p <- pt(abs(t), df = n - (1 + k), lower.tail = FALSE) * 2
    if (p == 0) {
      pvalue <- 1e-40
    } else {
      p <- p
    }
    pvalue[i.M] <- p
    beta[i.M] <- b
    sd[i.M] <- se
    Z[i.M] <- b/se
  }
  out <- data.frame(beta = beta, sd = sd, pvalue = pvalue, Z = Z)
  return(out)
}

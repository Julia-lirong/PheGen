

#' Simulate the matrix of genetic effects or noise effects
#'
#' @rdname GeneticNoiseEffect
#' @aliases geneticEffect noiseEffect
#'
#' @param N The number of samples to simulate.
#' @param M The number of variants in a region to simulate.
#' @param K The number of phenotypes to simulate.
#' @param W The weight assigned to each variant, default is 1.
#' @param V.g The covariance matrix among columns, default is 0.5 for
#'   non-diagnal element.
#' @param V.e The covariance matrix among rows, default is 0.5 for
#'   each non-diagnal element.
#' @param method String specifying the matrix decomposition used to determine
#'   the matrix root of sigma. Possible methods are eigenvalue decomposition
#'   ("eigen", default), singular value decomposition ("svd"), and
#'   Cholesky decomposition ("chol"). The Cholesky is typically fastest,
#'   not by much though.
#' @return \code{geneticEffect():} Matrix of genetic effects, rows correspond
#'    to variants, columns correspond to phenotypes. \code{noiseEffect():}
#'    Matrix of noise term, rows correspond to sample, columns correspond to
#'    phenotypes.
#' @importFrom mvtnorm rmvnorm
#' @export
#'
#' @examples
#' out <- geneticEffect(20, 10)


geneticEffect <- function(M, K, W = NULL, V.g = NULL, method = "eigen") {

  if (is.null(W)) {
    W <- matrix(0, nrow = M, ncol = M)
    diag(W) <- 1
  }
  if (dim(W)[1] != M) {
    stop("The row or column dimension of matrix W should equal to the
         number of varaints")
  }
  if (is.null(V.g)) {
    V.g <- matrix(0.5, nrow = K, ncol = K)
  }
  if (dim(V.g)[1] != K) {
    stop("The row or column dimension of matrix V.g should equal to the
         number of phenotypes")
  }
  # Check if the sigma is symmetric
  if (!isSymmetric(V.g)) {
    V.g = 0.5 * (V.g + t(V.g))
  }
  # Expand the matrix with each block being a diagonal matrix
  V.g <- kronecker(V.g, W)
  beta <- suppressWarnings(mvtnorm::rmvnorm(n = 1, sigma = V.g, method = method))
  beta <- matrix(beta, ncol = K)
  return(beta)
}


#' @rdname GeneticNoiseEffect
#' @export
#' @examples
#' out <- noiseEffect(20, 5)
noiseEffect <- function(N, K, V.e = NULL, method = "eigen") {

  if (is.null(V.e)) {
    V.e <- matrix(0.5, nrow = K, ncol = K)
  }
  if (dim(V.e)[1] != K) {
    stop("The row or column dimension of matrix V.e should equal to the
         number of phenotypes")
  }
  # Check if the sigma is symmetric
  if (!isSymmetric(V.e)) {
    V.e = 0.5 * (V.e + t(V.e))
  }
  beta <- suppressWarnings(mvtnorm::rmvnorm(n = N, sigma = V.e, method = method))
  return(beta)
}






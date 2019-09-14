#' @title Matrix bending to positive-definite
#'
#' @description Bending a non-positive-definite matrix to positive-definite, using weighted or unweighted methods.
#' 
#' @param inmat : The \code{matrix} to be bended.
#'
#' @param wtmat : The weight \code{matrix} for weighted bending. If no input is provided, the unweighted method (default) is used.
#'
#' @param reciprocal : If \code{TRUE}, reciprocal of the weighting factors are used. If no input is provided, default = \code{FALSE}.
#'
#' @param max.iter : Maximum number of iterations. If no input is provided, default = 10000.
#'
#' @param small.positive : A small positive value replacing smaller eigenvalues. If no input is provided, default = 0.0001.
#'
#' @return bended : The output bended \code{matrix}.
#'
#' @examples
#' # Test data
#' V = matrix(nrow=5, ncol=5, c( # matrix to be bended
#' 100,  95,  80,  40,  40,
#'  95, 100,  95,  80,  40,
#'  80,  95, 100,  95,  80,
#'  40,  80,  95, 100,  95,
#'  40,  40,  80,  95, 100))
#' W = matrix(nrow=5, ncol=5, c( # matrix of weights
#' 1000,  500,   20,   50,  200,
#'  500, 1000,  500,    5,   50,
#'   20,  500, 1000,   20,   20,
#'   50,    5,   20, 1000,  200,
#'  200,   50,   20,  200, 1000))
#'
#' # Example 1: Unweighted bending
#' bend(V)
#'
#' # Example 2: Weighted bending using reciprocal of the weighting factors
#' bend(inmat=V, wtmat=W, reciprocal=TRUE)
#'
#' # Example 3: Bending with fixed elements
#' ## Assume we want to keep V[1:2, 1:2] constant.
#' W2 = W; W2[1:2, 1:2] = 0
#' bend(inmat=V, wtmat=W2, reciprocal=TRUE)
#'
#' # Example 4: Bending a correlation matrix
#' V2 = cov2cor(V)
#' bend(V2, W, reciprocal=TRUE)
#'
#' @export
bend = function(inmat, wtmat, reciprocal=FALSE, max.iter=10000, small.positive=0.0001) {
  if(missing(wtmat)) wtmat = matrix(1, nrow=nrow(inmat), ncol=ncol(inmat))
  # Check for possible errors
  if(!is.matrix(inmat)) stop("ERROR: inmat object is not a matrix.")
  if(!isSymmetric(inmat)) stop("ERROR: The matrix is asymmetric.")
  if(!is.matrix(wtmat)) stop("ERROR: wtmat object is not a matrix.")
  if(!isSymmetric(wtmat)) stop("ERROR: The weight matrix is asymmetric.")
  if(!identical(ncol(inmat), ncol(wtmat))) stop("ERROR: The matrices are incompatible in size.")
  if(length(wtmat[wtmat < 0]) > 0) stop("ERROR: Found negative elements in the weight matrix.")
  if(!is.logical(reciprocal)) stop("ERROE: reciprocal should be a logical variable.")
  if(small.positive <= 0 | small.positive >= 0.1) stop("ERROR: 0 < small.positive < 0.1 is not met.")
  # Check if inmat is a correlation matrix
  if(identical(diag(inmat), rep(1, dim(inmat)[1]))) diag(wtmat) = 0
  # Get reciprocal of the weighting factors if needed
  if(reciprocal) {
    for(i in 1:nrow(wtmat))
    {
      for(j in 1:ncol(wtmat))
      {
        if(wtmat[i,j]!=0) wtmat[i,j] = 1/wtmat[i,j]
      }
    }
  }
  wtmat = wtmat/max(wtmat)
  # Start bending
  bended = inmat
  eigenval = eigen(bended, symmetric=TRUE)$values
  eigenvec = eigen(bended, symmetric=TRUE)$vectors
  if(length(eigenval[eigenval<=0])==0) stop("No action is needed. The matrix is positive-definite.")
  i = 0
  while(i < max.iter & length(eigenval[eigenval < 0]) > 0)
  {
    eigenval[eigenval < small.positive] <- small.positive
    UDtU = eigenvec %*% diag(eigenval) %*% t(eigenvec)
    bended = bended - (bended - UDtU) * wtmat
    eigenval = eigen(bended, symmetric=TRUE)$values
    eigenvec = eigen(bended, symmetric=TRUE)$vectors
    i = i + 1
    if(round(i %% 100)==0) message(paste('Iteration =', i))
  }
  if(length(eigenval[eigenval<=0]) > 0) {
    message(paste("WARNING: convergence was not met after", i, "iterations."))
  } else {
    message(paste("Convergence met after", i, "iterations."))
  }
  return(bended)
}

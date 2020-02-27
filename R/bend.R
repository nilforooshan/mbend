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
#' @param method : \code{"hj"} (Jorjani et al., 2003) or \code{"lrs"} (Schaeffer, 2010), default = \code{"hj"}
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
#' ## The default method (Jojani et al. 2003) is used.
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
#' # Example 5: Bending using the method of Schaeffer (2010)
#' bend(inmat=V, method="lrs")
#'
#' # Example 6: Bending a correlation matrix using the method of Schaeffer (2010)
#' bend(V2, method="lrs")
#'
#' # Example 7: Bending the same correlation matrix using a weighted development of Schaeffer (2010)
#' bend(V2, W, reciprocal=TRUE, method="lrs")
#'
#' # Example 8: Bending a covariance matrix using a weighted development of Schaeffer (2010)
#' bend(V, W, reciprocal=TRUE, method="lrs")
#'
#' @export
bend = function(inmat, wtmat, reciprocal=FALSE, max.iter=10000, small.positive=0.0001, method="hj") {
  if(!is.matrix(inmat)) stop("inmat object is not a matrix.")
  if(!identical(rownames(inmat), colnames(inmat))) stop("rownames and colnames of the input matrix are not identical.")
  if(!isSymmetric(inmat)) stop("The matrix is asymmetric.")
  bended = inmat
  eigenval = eigen(bended, symmetric=TRUE)$values
  # Report what you've got.
  if(missing(wtmat)) {
    message("Unweighted bending")
  } else {
    message("Weighted bending")
    message("reciprocal = ", reciprocal)
  }
  message("max.iter = ", max.iter)
  if(method=="hj") message("small.positive = ", small.positive)
  message("method = ", method)
  message("Initial eigenvalues:")
  for(i in eigenval) message("  ", i)
  if(length(eigenval[eigenval<=0]) > 0) { # If non-PD
    N = nrow(inmat)
    if(missing(wtmat)) {
      wtmat = matrix(1, nrow=N, ncol=N)
    } else {
      if(!is.matrix(wtmat)) stop("wtmat object is not a matrix.")
      if(!identical(rownames(wtmat), colnames(wtmat))) stop("rownames and colnames of the weight matrix are not identical.")
      if(!isSymmetric(wtmat)) stop("The weight matrix is asymmetric.")
      if(ncol(inmat)!=ncol(wtmat)) stop("The matrices are incompatible in size.")
      if(length(wtmat[wtmat <0]) > 0) stop("Found negative elements in the weight matrix.")
    }
    # Check if inmat is a correlation matrix
    correl = FALSE
    if(identical(diag(inmat), rep(1, N))) correl = TRUE
    if(correl) diag(wtmat) = 0
    if(!is.logical(reciprocal)) stop("reciprocal should be a logical variable.")
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
    if(!method %in% c("hj","lrs")) stop("Use method \"hj\" or \"lrs\".")
    if(method=="lrs") {
      message("Source: Schaeffer, L. R. (2010). http://animalbiosciences.uoguelph.ca/~lrs/piksLRS/PDforce.pdf")
    }
    eigenvec = eigen(bended, symmetric=TRUE)$vectors
    m = length(eigenval[eigenval < 0])
    i = 0
    while(i < max.iter & m > 0)
    {
      # Updating eigenvalues
      if(method=="hj") {
        if(small.positive <= 0) stop("small.positive is not positive.")
        if(i==0) {
          # If small.positive > smallest positive eigenvalue, replace it.
          sp.ev = eigenval[eigenval > 0]
          sp.ev = sp.ev[length(sp.ev)]
          if(small.positive > sp.ev) {
            small.positive = sp.ev/10
            message("small.positive was overwritten by ", small.positive)
          }
        }
        eigenval[eigenval < small.positive] <- small.positive
      } else { # method="lrs"
        v = as.numeric(eigenval < 0)
        S = sum(v*eigenval)*2
        W = (S*S*100)+1
        small.positive = eigenval[N - m]
        k = N - m + 1
        for(j in k:N) eigenval[j] = small.positive * (S-eigenval[j])*(S-eigenval[j])/W
        # Make sure that UDtU is PD.
        UDtU = eigenvec %*% diag(eigenval) %*% t(eigenvec)
        eigenval = eigen(UDtU, symmetric=TRUE)$values
        eigenval[eigenval < 0] <- small.positive/10
      }
      UDtU = eigenvec %*% diag(eigenval) %*% t(eigenvec)
      bended = bended - (bended - UDtU) * wtmat
      eigenval = eigen(bended, symmetric=TRUE)$values
      eigenvec = eigen(bended, symmetric=TRUE)$vectors
      m = length(eigenval[eigenval < 0])
      i = i + 1
      if(round(i %% 1000)==0) message(paste('Iteration =', i))
    }
    if(length(eigenval[eigenval < 0]) > 0) {
      message(paste("WARNING: convergence was not met after", i, "iterations."))
    } else {
      message(paste("Convergence met after", i, "iterations."))
    }
  } else {
    message("No action was required. The matrix is positive-definite.")
  }
  AD = abs(bended - inmat)
  minAD = which(AD==min(AD), arr.ind=TRUE)[1,]
  maxAD = which(AD==max(AD), arr.ind=TRUE)[1,]
  AD = AD[upper.tri(AD, diag=TRUE)]
  w = wtmat[upper.tri(wtmat, diag=TRUE)]
  cor.dat = data.frame(V1=inmat[upper.tri(inmat, diag=TRUE)], V2=bended[upper.tri(bended, diag=TRUE)])
  message("Final eigenvalues:")
  for(i in eigenval) message("  ", i)
  message("Minimum absolute deviation = ", min(AD), " [", minAD[1], ",", minAD[2], "]")
  message("Maximum absolute deviation = ", max(AD), " [", maxAD[1], ",", maxAD[2], "]")
  message("Upper.tri average absolute deviation = ", mean(AD))
  message("Upper.tri correlation = ", cor(cor.dat)[1,2])
  if(any(!w %in% 0:1)) {
    message(length(w), " Upper.tri elements; for the ", length(w[w>0]), " elements that are not set to constant:")
    message("  Upper.tri weighted average absolute deviation = ", sum(AD[w>0]/w[w>0])/sum(1/w[w>0]))
    message("  Upper.tri weighted correlation = ", cov.wt(cor.dat[w>0,], wt=1/w[w>0], cor=TRUE)$cor[1,2])
  }
  return(bended)
}


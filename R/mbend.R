#' @title Matrix Bending
#'
#' @docType package
#'
#' @name mbend-package
#'
#' @author Mohammad Ali Nilforooshan \email{m.a.nilforooshan@gmail.com}
#'
#' @description
#' (Co)variance or correlation matrices are required in multivariate mixed models. For example, in multi-trait animal models, genetic and residual (co)variance matrices are required, involving the traits of interest. These matrices need to be positive definite (PD) and invertible. Variance component estimation is computationally expensive, especially for big data and many variables. As a result, the full (co)variance matrix may be assembled by combining smaller (co)variance matrices from variance component estimation analyses on subsets of variables. Possible missing covariances are filled with values from the literature or the best possible guess. Consequently, the assembled matrix may not be PD, and it needs to be bent to a PD matrix before being used in the model.
#'
#' @details
#' A method for weighted bending of (co)variance matrices was developed by Jorjani et al. (2003), in which the matrix of interest is decomposed to matrices of eigenvectors and eigenvalues. Iteratively, eigenvalues smaller than a small possitive value (close to 0) are replaced with that small positive value and the matrix is rebuilt, until the convergence is met (i.e., all eigenvalues being positive). Because there are different amount of data and certainty associated with different elements of the matrix, wighting factors should be involved, which are introduced through a symmetric matrix. Certainty associated with the elements of the matrix and the corresponding weights are inversely related. For example, the reciprocal of the number of common data points (i.e., data points in common between pairs of variables) can be used as weighting factors. Alternatively, number of data points can be used directly by setting the argument \code{reciprocal = TRUE}.
#' To keep specific elements of the matrix unchanged during the bending process, set corresponding weights to 0. Providing no weight matrix is equivalent to unweighted bending.
#' Another method implemented in this package is from Schaeffer (2014), which can be defined by using the argument \code{method = "lrs"}. In this methos, negative eigenvalues are replaced with positive values in a descending order. If no method is defined, the default \code{method = "hj"} (Jorjani et al., 2003) is used. As a development to the method of Schaeffer (2014), a weight matrix can be used for weighted bending.
#' Any (co)variance matrix with all diagonal elements equal to 1 is considered as a correlation matrix by the program.
#'
#' @references
#' Jorjani, H., et al. (2003). A Simple Method for Weighted Bending of Genetic (Co)variance Matrices. J. Dairy Sci., 86:677-679. <doi:10.3168/jds.S0022-0302(03)73646-7>
#'
#' Schaeffer, L. R. (2014). Making covariance matrices positive definite.
#' Available at: \href{http://animalbiosciences.uoguelph.ca/~lrs/ELARES/PDforce.pdf}{Link}
NULL

# mbend: An R package for bending non-positive-definite matrices to positive-definite

Mohammad Ali Nilforooshan

<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0003-0339-5442" href="https://orcid.org/0000-0003-0339-5442" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0003-0339-5442</a></div>

*14 September 2019*

### Installation

```r
# install.packages("devtools")
devtools::install_github('nilforooshan/mbend')
```

Alternatively:

```r
path = "https://github.com/nilforooshan/Link-resources/raw/master/Resources/"
installer = file.path(tempdir(), "mbend.tar.gz")
download.file(paste0(path, "mbend_1.0.0.tar.gz"), destfile=installer)
install.packages(installer, repos=NULL, type='source')
```

[PDF manual](https://github.com/nilforooshan/mbend/blob/master/mbend-manual.pdf)

## Summary

In statistical multivariate mixed models, the (co)variance structure between different variables are required for different random effects. For example, in the context of animal breeding and genetics, genetic and residual (co)variances between different traits in the multi-trait analysis are needed. With additional random effects in the model, such as maternal genetic and permanent environmental effects, additional (co)variance matrices are needed. Variance component estimation for large number of traits and big data is computationally challenging. Therefore, a full (co)variance matrix for large number of traits is usually assembled from smaller (co)variance matrices obtained from several variance component estimation analyses, including different subsets of traits. Finally, possible missing elements are filled with values from the literature or the best possible guess. The resulting matrix is generally not positive-definite (PD), and because the inverse of the (co)variance matrices are required in the model, these matrices need to be PD.

An iterative procedure called "bending" is used to convert a non-PD matrix to a PD matrix closest from the non-PD matrix. In an "unweighted" bending procedure, there is no precision associated with matrix elements, and elements have equal probability of changing during the bending procedure (Jorjani et al., 2003), which is not preferable. Jorjani et al. (2003) introduced "weighted" bending. Because the precision of matrix elements and the corresponding weights are inversely related, they used the reciprocal of the number of data points (in common between pairs of traits) as weights.

``mbend`` is an R package for bending (weighted or unweighted) a non-PD matrix to PD.

## Application

Though, the application of the weighted bending method is given for animal breeding and genetics, and multivariate mixed models in general, it can be applicable in any field of science, where transformation of a non-PD matrix to a PD matrix is required.

For the examples, the same input matrices from Jorjani et al. (2003) are used, where `V` is the matrix to be bended, and `W` is the matrix of weights. Bending `V` with reciprocals of `W`:

```r
V = matrix(nrow=5, ncol=5, c( # matrix to be bended
  100,  95,  80,  40,  40,
   95, 100,  95,  80,  40,
   80,  95, 100,  95,  80,
   40,  80,  95, 100,  95,
   40,  40,  80,  95, 100))
W = matrix(nrow=5, ncol=5, c( # matrix of weights
  1000,  500,   20,   50,  200,
   500, 1000,  500,    5,   50,
    20,  500, 1000,   20,   20,
    50,    5,   20, 1000,  200,
   200,   50,   20,  200, 1000))
bend(inmat=V, wtmat=W, reciprocal=TRUE)
# Iteration = 100
# Iteration = 200
# Iteration = 300
# Iteration = 400
# Convergence met after 428 iterations.
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 100.16154  94.51747  82.93843  43.56681  39.18292
# [2,]  94.51747 100.62488  93.98150  59.99921  45.84555
# [3,]  82.93843  93.98150 100.69547  84.89294  73.13431
# [4,]  43.56681  59.99921  84.89294 100.30697  94.23170
# [5,]  39.18292  45.84555  73.13431  94.23170 100.17942
```

Please note that the default parameters `max.iter = 10000, small.positive = 0.0001` could be omitted. Comparing the results with unweighted bending:

```r
bend(V)
# Convergence met after 1 iterations.
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 103.16918  90.83313  79.47122  44.53731  37.07254
# [2,]  90.83313 106.49734  94.18961  74.07039  44.53731
# [3,]  79.47122  94.18961 102.31355  94.18961  79.47122
# [4,]  44.53731  74.07039  94.18961 106.49734  90.83313
# [5,]  37.07254  44.53731  79.47122  90.83313 103.16918
```

In comparison with Jorjani et al. (2003), weights are divided by the maximum weight, which makes the convergence faster. The program considers any (co)variance matrix with all diagonal values equal to 1, as a correlation matrix. It takes a different approach to bend correlation matrices compared to Jorjani et al. (2003), by forcing the diagonal values to 1 during the bending procedure. For example:

```r
V2 = cov2cor(V)
bend(V2, W, reciprocal=TRUE)
# Iteration = 100
# Iteration = 200
# Convergence met after 286 iterations.
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 1.0000000 0.9447536 0.8342673 0.4377103 0.3912583
# [2,] 0.9447536 1.0000000 0.9385087 0.6004642 0.4630373
# [3,] 0.8342673 0.9385087 1.0000000 0.8394178 0.7248761
# [4,] 0.4377103 0.6004642 0.8394178 1.0000000 0.9418815
# [5,] 0.3912583 0.4630373 0.7248761 0.9418815 1.0000000
```

In order to force bending without changing specific elements of the matrix, the corresponding weights must be set to 0 (i.e., full precision). For example, to keep the first 2 &times; 2 block of `V` unchanged during the bending procedure:

```r
W2 = W; W2[1:2, 1:2] = 0
bend(V, W2, reciprocal=TRUE)
# Iteration = 100
# Iteration = 200
# Iteration = 300
# Iteration = 400
# Convergence met after 475 iterations.
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 100.00000  95.00000  83.70895  43.70489  39.12852
# [2,]  95.00000 100.00000  93.89940  59.58422  46.15038
# [3,]  83.70895  93.89940 100.73015  84.47456  72.91450
# [4,]  43.70489  59.58422  84.47456 100.31566  94.21299
# [5,]  39.12852  46.15038  72.91450  94.21299 100.18322
```

## References

Jorjani, H., Klie. L., & Emanuelson, U. (2000). A simple method for weighted bending of genetic (co)variance matrices. *J. Dairy Sci. 86(2)*: 677–679. doi:[10.3168/jds.S0022-0302(03)73646-7](https://doi.org/10.3168/jds.S0022-0302(03)73646-7)

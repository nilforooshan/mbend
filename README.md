# mbend: An R package for bending non-positive-definite matrices to positive-definite

Mohammad Ali Nilforooshan

<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0003-0339-5442" href="https://orcid.org/0000-0003-0339-5442" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0003-0339-5442</a></div>

---

## Description

*Bending non-positive-definite (symmetric) matrices to positive-definite, using weighted and unweighted methods*

The `mbend` package is used for bending a non-positive-definite matrix to a positive-definite (PD) matrix. For a matrix to be invertible, it has to be PD. The methods of Jorjani et al. (2003) and Schaeffer (2010) are used in this package. The unweighted method of Schaeffer (2010) is also extended to weighted bending, with the possibility of choosing between unweighted and weighted bending.

## Application

Start with loading the package library:

```r
library(mbend)
```

Consider the following non-PD matrix (Jorjani et al., 2003):

```r
V = matrix(nrow=5, ncol=5, c(
  100,  95,  80,  40,  40,
   95, 100,  95,  80,  40,
   80,  95, 100,  95,  80,
   40,  80,  95, 100,  95,
   40,  40,  80,  95, 100))
```

The following command can be used to bend matrix `V` to a PD matrix:

```r
bend(V)
#> Convergence met after 1 iterations.
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 103.16918  90.83313  79.47122  44.53731  37.07254
#> [2,]  90.83313 106.49734  94.18961  74.07039  44.53731
#> [3,]  79.47122  94.18961 102.31355  94.18961  79.47122
#> [4,]  44.53731  74.07039  94.18961 106.49734  90.83313
#> [5,]  37.07254  44.53731  79.47122  90.83313 103.16918
```

The above command is equivalent to `bend(inmat=V)`, where `inmat` is the argument that takes the matrix to be bended, or the following command:

```r # eval=FALSE
bend(V, max.iter=10000, small.positive=0.0001, method="hj")
```

This runs the unweighted bending method of Jorjani et al. (2003) (`method="hj"`), with maximum 10000 number of iterations (`max.iter=10000`), and eigenvalues smaller than 0.0001 are replaced with this small positive value (`small.positive=0.0001`). Providing the default parameters, the corresponding arguments can be omitted (e.g., `bend(V)`).

There might be different precision involved with different elements of a non-PD matrix. In this case, a weighted bending is recommended. Jorjani et al. (2003) used the reciprocal of the number data points in common between pairs of variables, as weights. Considering the following matrix for the number of data points in common between variables (Jorjani et al., 2003):

```r
W = matrix(nrow=5, ncol=5, c(
  1000,  500,   20,   50,  200,
   500, 1000,  500,    5,   50,
    20,  500, 1000,   20,   20,
    50,    5,   20, 1000,  200,
   200,   50,   20,  200, 1000))
```

Matrix `V` is bended using the following command:

```r
bend(inmat=V, wtmat=W, reciprocal=TRUE)
#> Iteration = 100
#> Iteration = 200
#> Iteration = 300
#> Iteration = 400
#> Convergence met after 428 iterations.
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 100.16154  94.51747  82.93843  43.56681  39.18292
#> [2,]  94.51747 100.62488  93.98150  59.99921  45.84555
#> [3,]  82.93843  93.98150 100.69547  84.89294  73.13431
#> [4,]  43.56681  59.99921  84.89294 100.30697  94.23170
#> [5,]  39.18292  45.84555  73.13431  94.23170 100.17942
```

Using `wtmat=1/W, reciprocal=FALSE`, the argument `reciprocal` could be omitted, because `FALSE` is the default parameter for the argument `reciprocal`. For the same reason, `max.iter=10000, small.positive=0.0001, method="hj"` are omitted, unless different parameters are provided to these arguments.

If there is high confidence about some elements of the non-PD matrix to remain unchanged after bending, the corresponding weights are set to zero. For example, to keep the first 2 &times; 2 block of `V` unchanged during the bending procedure:

```r
W2 = W; W2[1:2, 1:2] = 0
bend(V, W2, reciprocal=TRUE)
#> Iteration = 100
#> Iteration = 200
#> Iteration = 300
#> Iteration = 400
#> Convergence met after 475 iterations.
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 100.00000  95.00000  83.70895  43.70489  39.12852
#> [2,]  95.00000 100.00000  93.89940  59.58422  46.15038
#> [3,]  83.70895  93.89940 100.73015  84.47456  72.91450
#> [4,]  43.70489  59.58422  84.47456 100.31566  94.21299
#> [5,]  39.12852  46.15038  72.91450  94.21299 100.18322
```

To bend `V` using the method of Schaeffer (2010):

```r
bend(inmat=V, method="lrs")
#> Source: Schaeffer, L. R. (2010). http://animalbiosciences.uoguelph.ca/~lrs/piksLRS/PDforce.pdf
#> NOTE: argument small.positive is overwritten.
#> Convergence met after 1 iterations.
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 103.18978  90.82704  79.43676  44.56754  37.06769
#> [2,]  90.82704 106.54177  94.13680  74.06295  44.56754
#> [3,]  79.43676  94.13680 102.46429  94.13680  79.43676
#> [4,]  44.56754  74.06295  94.13680 106.54177  90.82704
#> [5,]  37.06769  44.56754  79.43676  90.82704 103.18978
```

The method of Schaeffer (2010) does not require the argument `small.positive`, and this argument is ignored. This method is originally an unweighted bending method. However, in this package, it is extended to accommodate weighted bending (i.e., a combination of Schaeffer (2010) and Jorjani et al. (2003) methods). Weighted bending of `V` with reciprocals of `W` using the method of Schaeffer (2010):

```r
bend(V, W, reciprocal=TRUE, method="lrs")
#> Source: Schaeffer, L. R. (2010). http://animalbiosciences.uoguelph.ca/~lrs/piksLRS/PDforce.pdf
#> NOTE: argument small.positive is overwritten.
#> Iteration = 100
#> Iteration = 200
#> Iteration = 300
#> Iteration = 400
#> Iteration = 500
#> Iteration = 600
#> Iteration = 700
#> Convergence met after 787 iterations.
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 100.16194  94.51554  82.97457  43.57119  39.18071
#> [2,]  94.51554 100.62754  93.97669  60.00766  45.86380
#> [3,]  82.97457  93.97669 100.69806  84.86261  73.10529
#> [4,]  43.57119  60.00766  84.86261 100.30785  94.22946
#> [5,]  39.18071  45.86380  73.10529  94.22946 100.18003
```

Function `bend` automatically considers any matrix with all diagonal elements equal to one, as a correlation matrix. Considering the correlation matrix `V2` (`V` converted to a correlation matrix):

```r
V2 = cov2cor(V)
bend(V2, W, reciprocal=TRUE)
#> Iteration = 100
#> Iteration = 200
#> Convergence met after 286 iterations.
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.9447536 0.8342673 0.4377103 0.3912583
#> [2,] 0.9447536 1.0000000 0.9385087 0.6004642 0.4630373
#> [3,] 0.8342673 0.9385087 1.0000000 0.8394178 0.7248761
#> [4,] 0.4377103 0.6004642 0.8394178 1.0000000 0.9418815
#> [5,] 0.3912583 0.4630373 0.7248761 0.9418815 1.0000000
```

Because the argument `method` is not provided, the default parameter `"hj"` is used. To do the same using the method of Schaeffer (2010):

```r
bend(V2, W, reciprocal=TRUE, method="lrs")
#> Source: Schaeffer, L. R. (2010). http://animalbiosciences.uoguelph.ca/~lrs/piksLRS/PDforce.pdf
#> NOTE: argument small.positive is overwritten.
#> Iteration = 100
#> Iteration = 200
#> Iteration = 300
#> Iteration = 400
#> Iteration = 500
#> Iteration = 600
#> Iteration = 700
#> Convergence met after 700 iterations.
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.9447383 0.8345316 0.4377385 0.3912447
#> [2,] 0.9447383 1.0000000 0.9384618 0.6005923 0.4631893
#> [3,] 0.8345316 0.9384618 1.0000000 0.8392102 0.7245638
#> [4,] 0.4377385 0.6005923 0.8392102 1.0000000 0.9418734
#> [5,] 0.3912447 0.4631893 0.7245638 0.9418734 1.0000000
```

Bending the same correlation matrix using the unweighted Schaeffer (2010):

```r
bend(V2, method="lrs")
#> Source: Schaeffer, L. R. (2010). http://animalbiosciences.uoguelph.ca/~lrs/piksLRS/PDforce.pdf
#> NOTE: argument small.positive is overwritten.
#> Convergence met after 39 iterations.
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.8935779 0.7794136 0.4716892 0.3599505
#> [2,] 0.8935779 1.0000000 0.9128285 0.7211811 0.4716892
#> [3,] 0.7794136 0.9128285 1.0000000 0.9128285 0.7794136
#> [4,] 0.4716892 0.7211811 0.9128285 1.0000000 0.8935779
#> [5,] 0.3599505 0.4716892 0.7794136 0.8935779 1.0000000
```

---

## References

Jorjani, H., Klie. L., & Emanuelson, U. (2000). A simple method for weighted bending of genetic (co)variance matrices. *J. Dairy Sci.* 86(2): 677--679. [doi:10.3168/jds.S0022-0302(03)73646-7](https://doi.org/10.3168/jds.S0022-0302(03)73646-7)

Schaeffer, L. R. (2010). Modification of negative eigenvalues to create positive definite matrices and approximation of standard errors of correlation estimates. Available at: http://animalbiosciences.uoguelph.ca/~lrs/piksLRS/PDforce.pdf

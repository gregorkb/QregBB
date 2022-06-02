
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QregBB

<!-- badges: start -->

<!-- badges: end -->

The R package QregBB accompanies the paper:

Gregory, K.B., Lahiri, S.N., Nordman, D.J. (2018). A smooth block
bootstrap for quantile regression with time series. *Annals of
Statistics* 46(3), 1138-1166

Find the paper at <https://projecteuclid.org/euclid.aos/1525313078>.

## Installation

You can install the development version of QregBB from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gregorkb/QregBB")
```

## Example

The main function in the package is `QregBB`, which performs the MBB,
SMBB, ETBB, and SETBB bootstrap procedures (all at once) for estimating
the sampling distributions of quantile regression estimators with time
series data.

``` r
library(QregBB)

n <- 50
X1 <- arima.sim(model=list(ar=c(.7,.1)),n)
X2 <- arima.sim(model=list(ar=c(.2,.1)),n)
e <- arima.sim(model=list(ar=c(.7,.1)),n)
Y <- X1 + e
X <- cbind(rep(1,n),X1,X2)
QregBB.out <- QregBB(Y,X,tau=.5,l=4,B=500,h=NULL,alpha=0.05)

QregBB.out
#> 
#> Call:
#> QregBB(Y = Y, X = X, tau = 0.5, l = 4, B = 500, h = NULL, alpha = 0.05)
#> 
#> Coefficients:
#>             Estimate SE (MBB) SE (SMBB) SE (ETBB) SE (SETBB)
#> (Intercept)  1.01493  0.68347   0.62171   0.62154    0.59565
#> beta_1       0.79606  0.32682   0.29734   0.30377    0.28419
#> beta_2       0.31152  0.34370   0.38398   0.39091    0.37435
#> 
#> Confidence intervals:
#>             Estimate lower (MBB) upper (MBB) lower (SMBB) upper (SMBB)
#> (Intercept)  1.01493    -0.27172     2.30680     -0.14704      2.29554
#> beta_1       0.79606     0.00708     1.07686      0.16048      1.35850
#> beta_2       0.31152    -0.38721     0.91367     -0.45135      1.01176
#>             lower (ETBB) upper (ETBB) lower (SETBB) upper (SETBB)
#> (Intercept)     -0.10149      2.30319      -0.03575       2.35307
#> beta_1           0.14920      1.16432       0.16888       1.31251
#> beta_2          -0.72947      0.66229      -0.48280       0.97010
```

The function `getNPPIblksizesQR` implements the block size selection
method described in Gregory et al. (2018) for MBB, SMBB, ETBB, and
SETBB.

``` r
blksize.out <- getNPPIblksizesQR(Y,X,tau=.5)
blksize.out
#> $l.opt.MBB
#> [1] 25
#> 
#> $l.opt.ETBB
#> [1] 6
#> 
#> $l.opt.SMBB
#> [1] 1
#> 
#> $l.opt.SETBB
#> [1] 2
```

<!--
You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
-->

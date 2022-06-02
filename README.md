
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

## Examples

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
#> (Intercept) -1.19048  0.30832   0.35614   0.29801    0.34243
#> beta_1       1.37104  0.14455   0.16978   0.13906    0.16874
#> beta_2       0.00682  0.18293   0.24371   0.20035    0.24457
#> 
#> Confidence intervals:
#>             Estimate lower (MBB) upper (MBB) lower (SMBB) upper (SMBB)
#> (Intercept) -1.19048    -2.01656    -0.78747     -1.93770     -0.61313
#> beta_1       1.37104     1.00985     1.56780      1.00313      1.64540
#> beta_2       0.00682    -0.44854     0.28923     -0.47015      0.45646
#>             lower (ETBB) upper (ETBB) lower (SETBB) upper (SETBB)
#> (Intercept)     -1.97114     -0.76846      -1.83370      -0.58475
#> beta_1           1.03540      1.59665       0.96382       1.67961
#> beta_2          -0.49587      0.29042      -0.47396       0.43994
```

The function `getNPPIblksizesQR` implements the block size selection
method described in Gregory et al. (2018) for MBB, SMBB, ETBB, and
SETBB.

``` r
blksize.out <- getNPPIblksizesQR(Y,X,tau=.5)
blksize.out
#> $l.opt.MBB
#> [1] 23
#> 
#> $l.opt.ETBB
#> [1] 3
#> 
#> $l.opt.SMBB
#> [1] 25
#> 
#> $l.opt.SETBB
#> [1] 3
```

<!--
You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>.

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
-->

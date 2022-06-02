\name{QregBB-package}
\alias{QregBB-package}
\docType{package}
\title{
\packageTitle{QregBB}
}
\description{
\packageDescription{QregBB}
}
\details{
Implements moving-blocks bootstrap and extended tapered-blocks bootstrap, as well as smooth versions of each, for quantile regression in time series. This package accompanies Gregory et al. (2018).
\packageDESCRIPTION{QregBB}
\packageIndices{QregBB}
The main function is the \code{QregBB} function, which implements the moving-blocks bootstrap (MBB), the extended tapered-blocks bootstrap (ETBB), and smooth versions of each (SMBB, SETBB). The function \code{getNPPIblksizesQR} chooses the block size based on the non-parametric plug-in method described in Lahiri (2013). For the smooth methods, the bandwidth is chosen by using the function \code{bw.SJ} function on the fitted residuals; then the bandwidth matrix is the identity matrix times the value returned by \code{bw.SJ}.
}
\author{
\packageAuthor{QregBB}

Maintainer: \packageMaintainer{QregBB}
}
\references{

Gregory, K. B., Lahiri, S. N., & Nordman, D. J. (2018). A smooth block bootstrap for quantile regression with time series. \emph{The Annals of Statistics}, 46(3), 1138-1166.

Lahiri, S. N. (2013). Resampling methods for dependent data. Springer Science & Business Media.

}


\keyword{ block bootstrap }
\examples{
n <- 100
X1 <- arima.sim(model=list(ar=c(.7,.1)),n)
X2 <- arima.sim(model=list(ar=c(.2,.1)),n)
e <- arima.sim(model=list(ar=c(.7,.1)),n)
Y <- X1 + e
X <- cbind(rep(1,n),X1,X2)

QregBB.out <- QregBB(Y,X,tau=.5,l=4,B=500,h=NULL,alpha=0.05)
QregBB.out

}

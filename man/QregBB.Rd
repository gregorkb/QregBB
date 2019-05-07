\name{QregBB}
\alias{QregBB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QregBB(y = Y, x = X, tau, l = 5, B = 500, h = NULL, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{l}{
%%     ~~Describe \code{l} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y = Y, x = X, tau, l = 5, B = 500, h = NULL, alpha = 0.05) 
{
    n <- length(Y)
    d <- ncol(X) - 1
    model <- rq(Y ~ X[, -1], tau = tau)
    beta.hat <- model$coef
    if (length(h) == 0) {
        h <- bw.SJ(model$resid)
    }
    AllBBgetweights.out <- AllBBgetweights(n, l, B)
    MBB.weights <- AllBBgetweights.out$MBBweights
    ETBB.weights <- AllBBgetweights.out$ETBBweights
    m.l <- AllBBgetweights.out$m_l
    pi.tilde.MBB <- get.pi.tilde.MBB(n, l)
    pi.tilde.ETBB <- get.pi.tilde.ETBB(n, l)
    MBB.beta.tilde <- rq(Y ~ X[, -1], tau = tau, weights = pi.tilde.MBB)$coef
    ETBB.beta.tilde <- rq(Y ~ X[, -1], tau = tau, weights = pi.tilde.ETBB)$coef
    SMBB.beta.tilde <- nlm(f = S.tilde, p = beta.hat, Y = Y, 
        X = X, pi.tilde = pi.tilde.MBB, h = h, tau = tau, iterlim = 200)$estimate
    SETBB.beta.tilde <- nlm(f = S.tilde, p = beta.hat, Y = Y, 
        X = X, pi.tilde = pi.tilde.ETBB, h = h, tau = tau, iterlim = 200)$estimate
    MBB.beta.hat <- ETBB.beta.hat <- SMBB.beta.hat <- SETBB.beta.hat <- matrix(NA, 
        B, d + 1)
    for (b in 1:B) {
        MBB.beta.hat[b, ] <- rq(Y ~ X[, -1], tau = tau, weights = MBB.weights[, 
            b])$coef
        ETBB.beta.hat[b, ] <- rq(Y ~ X[, -1], tau = tau, weights = ETBB.weights[, 
            b])$coef
        Y.smooth <- Y + rnorm(n, 0, h)
        X.smooth <- X + cbind(rep(0, n), matrix(rnorm(n * d, 
            0, h), n, d))
        SMBB.beta.hat[b, ] <- rq(Y.smooth ~ X.smooth[, -1], tau = tau, 
            weights = MBB.weights[, b])$coef
        SETBB.beta.hat[b, ] <- rq(Y.smooth ~ X.smooth[, -1], 
            tau = tau, weights = ETBB.weights[, b])$coef
    }
    MBB.cov.est <- cov(MBB.beta.hat)
    ETBB.cov.est <- m.l * cov(ETBB.beta.hat)
    SMBB.cov.est <- cov(SMBB.beta.hat)
    SETBB.cov.est <- m.l * cov(SETBB.beta.hat)
    MBB.pivot <- sqrt(n) * (MBB.beta.hat - matrix(MBB.beta.tilde, 
        B, d + 1, byrow = TRUE))
    ETBB.pivot <- sqrt(m.l) * sqrt(n) * (ETBB.beta.hat - matrix(ETBB.beta.tilde, 
        B, d + 1, byrow = TRUE))
    SMBB.pivot <- sqrt(n) * (SMBB.beta.hat - matrix(SMBB.beta.tilde, 
        B, d + 1, byrow = TRUE))
    SETBB.pivot <- sqrt(m.l) * sqrt(n) * (SETBB.beta.hat - matrix(SETBB.beta.tilde, 
        B, d + 1, byrow = TRUE))
    MBB.loci <- beta.hat - apply(MBB.pivot, 2, quantile, 1 - 
        alpha/2)/sqrt(n)
    MBB.upci <- beta.hat - apply(MBB.pivot, 2, quantile, alpha/2)/sqrt(n)
    MBB.confint <- cbind(MBB.loci, MBB.upci)
    ETBB.loci <- beta.hat - apply(ETBB.pivot, 2, quantile, 1 - 
        alpha/2)/sqrt(n)
    ETBB.upci <- beta.hat - apply(ETBB.pivot, 2, quantile, alpha/2)/sqrt(n)
    ETBB.confint <- cbind(ETBB.loci, ETBB.upci)
    SMBB.loci <- beta.hat - apply(SMBB.pivot, 2, quantile, 1 - 
        alpha/2)/sqrt(n)
    SMBB.upci <- beta.hat - apply(SMBB.pivot, 2, quantile, alpha/2)/sqrt(n)
    SMBB.confint <- cbind(SMBB.loci, SMBB.upci)
    SETBB.loci <- beta.hat - apply(SETBB.pivot, 2, quantile, 
        1 - alpha/2)/sqrt(n)
    SETBB.upci <- beta.hat - apply(SETBB.pivot, 2, quantile, 
        alpha/2)/sqrt(n)
    SETBB.confint <- cbind(SETBB.loci, SETBB.upci)
    output <- list(MBB.pivot = MBB.pivot, ETBB.pivot = MBB.pivot, 
        SMBB.pivot = SMBB.pivot, SETBB.pivot = SETBB.pivot, MBB.confint = MBB.confint, 
        ETBB.confint = ETBB.confint, SMBB.confint = SMBB.confint, 
        SETBB.confint = SETBB.confint, MBB.cov.est = MBB.cov.est, 
        ETBB.cov.est = ETBB.cov.est, SMBB.cov.est = SMBB.cov.est, 
        SETBB.cov.est = SETBB.cov.est)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

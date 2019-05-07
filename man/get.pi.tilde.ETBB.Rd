\name{get.pi.tilde.ETBB}
\alias{get.pi.tilde.ETBB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.pi.tilde.ETBB(n, l, c = 0.43)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{l}{
%%     ~~Describe \code{l} here~~
}
  \item{c}{
%%     ~~Describe \code{c} here~~
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
function (n, l, c = 0.43) 
{
    if (l == 1) {
        return(rep(1/n, n))
    }
    else {
        omega.1 <- sum(omega((1:l - 0.5)/l, c = c))
        pi.tilde <- numeric(n)
        for (t in 1:(l - 1)) {
            pi.tilde[t] <- pi.tilde[n - t + 1] <- sum(omega((1:t - 
                0.5)/l, c = c)/omega.1)/(n - l + 1)
        }
        pi.tilde[l:(n - l + 1)] <- 1/(n - l + 1)
        return(pi.tilde)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

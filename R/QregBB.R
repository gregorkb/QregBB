

#' Retrieve weights from from Monte-Carlo block samples for the MBB and ETBB
#'
#' @param n length of observed time series.
#' @param l block length.
#' @param B the number of Monte-Carlo draws.
#' @param c a constant controlling the tapering of weight assignments over the blocks.
#' @param seed a random number generating seed (optional).
#' @return a list containing \code{n} by \code{B} matrices containing the the mass placed on the \code{n} observations by the MBB and ETBB empirical distributions in each of the \code{B} Monte-Carlo draws, the sampled indices for the block starting points, and the scalar \eqn{m_l} corresponding to the weight function and the block length.
#' @noRd
BBgetweights <- function(n,l,B,seed=NA,c=.43){

	if(is.na(seed)) seed <- floor(runif(1)*100000)

	BBweights <- .C("BBgetweights",
	n = as.integer(n),
	l = as.integer(l),
	B = as.integer(B),
	blkpts = integer(B*floor(n/l)),
	f_MBB = double(B*n+l),
	f_ETBB = double(B*n+l),
	as.integer(seed),
	weights = double(l),
	as.double(c),
	m_l = double(1))		

	ETBBweights <- matrix( BBweights$f_ETBB[1:(n*B)],n,B)
	MBBweights <- matrix( BBweights$f_MBB[1:(n*B)],n,B)
	blkpoints <-  matrix( BBweights$blkpts,floor(n/l),B) + 1

	output <-  list(ETBBweights = ETBBweights,
	                MBBweights = MBBweights,
	                blkpoints = blkpoints,
	                m_l = BBweights$m_l )
	
	return(output)
	
}


#' Function defining the tapered assignment of weights across the ETBB blocks
#'
#' @param u a real number.
#' @param c a constant controlling the tapering.
#' @return the value of the tapering function.
#' @noRd
omega <- function( u,c=.43){ ( u <= c) * u/c + ((u > c) & (u < 1 - c)) + (u >= 1-c) *(1 - u)/c}

#' Computes expected value of the bootstrap weights assigned to each time point by the ETBB
#'
#' @param n the length of the observed time series.
#' @param l the block length.
#' @param c a number controlling the tapering.
#' @return a vector of length \code{n} containing the expected values of the weights assigned to each time point by the ETBB.
#' @noRd
get.pi.tilde.ETBB <- function(n,l,c=.43)
{
	if(l==1)
	{
		return(rep(1/n,n))
		
	} else {
	
		omega.vec <- omega( (1:l-.5) / l , c = c)
		omega.sum <- sum( omega.vec )
		
		pi.tilde <- numeric(n)
		for(t in 1:(l-1))
		{
			pi.tilde[t] <- pi.tilde[n-t+1] <- sum(omega.vec / omega.sum ) / (n - l + 1)
		}
		
		pi.tilde[l:(n-l+1)] <- 1/(n - l + 1)
			
		return(pi.tilde)
	}
}


#' Computes expected value of the bootstrap weights assigned to each time point by the MBB
#'
#' @param n the length of the observed time series.
#' @param l the block length.
#' @return a vector of length \code{n} containing the expected values of the weights assigned to each time point by the MBB.
#' @noRd
get.pi.tilde.MBB <- function(n,l)
{
	if(l==1)
	{
		return(rep(1/n,n))
		
	} else {
			
		pi.tilde <- numeric(n)
		for(t in 1:(l-1))
		{
			pi.tilde[t] <- pi.tilde[n-t+1] <- (t/l)/(n - l + 1)
		}
		
		pi.tilde[l:(n-l+1)] <- 1/(n - l + 1)
			
		return(pi.tilde)
	
	}
}


#' Computes bootstrap expected value of the quantile regression objective function
#'
#' @param beta parameter vector at which to compute the objective function.
#' @param Y the vector of response values.
#' @param X the design matrix (including a column of ones for the intercept).
#' @param pi.tilde weights retrieved from the \code{get.pi.tilde.MBB} of \code{get.pi.tilde.ETBB} functions.
#' @param h a scalar bandwidth (bandwidth matrix is \code{h} times identity).
#' @param tau the quantile of interest.
#' @return the value of the bootstrap expected value of the quantile regression objective function at the given beta.
#' This function will be numerically maximized in order to define the SMBB and SETBB versions of the true parameter vector; that is, the maximizer of this function is used as the centering \eqn{\tilde \beta_n} in the pivot \eqn{\sqrt{n}(\hat \beta^*_n - \tilde \beta_n)} rather than the sample estimator.
#' @noRd
S.tilde <- function(beta,Y,X,pi.tilde,h,tau)
{
	
	var.eta <- h^2*(1 + sum(beta[-1]^2))
	e.hat <- Y - X %*% beta
	Phi.e.hat <- pnorm(e.hat,0,sqrt(var.eta))
	phi.e.hat <- dnorm(e.hat,0,sd=sqrt(var.eta))
	S.tilde <- sum(pi.tilde*( e.hat*(tau-(1-Phi.e.hat)) + var.eta*phi.e.hat) )
	
	return(S.tilde)

}

#' Implements MBB, ETBB, SMBB, and SETBB for quantile regression
#'
#' @param Y the vector of response values.
#' @param X the design matrix (including a column of ones for the intercept).
#' @param tau the quantile of interest.
#' @param l block size.
#' @param B the number of Monte Carlo bootstrap samples to draw.
#' @param h a scalar bandwidth (bandwidth matrix is \code{h} times identity).
#' @param alpha a significance level to which the returned confidence intervals will correspond.
#' @return A list is returned containing for the MBB, SMBB, ETBB, and SETBB the set of Monte Carlo draws of the pivot quantity \eqn{\sqrt{n}(\hat \beta^*_n - \tilde \beta_n)}, confidence intervals for each component of \eqn{\beta} corresponding to the specified confidence level, and estimates of the asymptotic covariance matrix of the pivot quantity \eqn{\sqrt{n}(\hat \beta_n - \beta)}.
#'
#' @seealso  A `print.QregBB` method exists which prints to the console the bootstrap standard errors for each coefficient estimator from the MBB, SMBB, ETBB, and SETBB methods as well as confidence intervals for each coefficient at the specified level.
#'
#' @references 
#' 
#' #' @references 
#' 
#' Gregory, K. B., Lahiri, S. N., & Nordman, D. J. (2018). A smooth block bootstrap for quantile regression with time series. \emph{The Annals of Statistics}, 46(3), 1138-1166.
#' 
#' @examples
#' # generate some data and perform block-bootstrap methods
#' n <- 100
#' X1 <- arima.sim(model=list(ar=c(.7,.1)),n)
#' X2 <- arima.sim(model=list(ar=c(.2,.1)),n)
#' e <- arima.sim(model=list(ar=c(.7,.1)),n)
#' Y <- X1 + e
#' X <- cbind(rep(1,n),X1,X2)
#'
#' QregBB.out <- QregBB(Y,X,tau=.5,l=4,B=500,h=NULL,alpha=0.05)
#' QregBB.out
#' @export
QregBB <- function(Y,X,tau,l,B=500,h=NULL,alpha=0.05)
{
	
	n <- length(Y)
	d <- ncol(X) - 1

	model <- rq(Y ~ X[,-1], tau = tau)
	beta.hat <- model$coef
	names(beta.hat) <- c("(Intercept)",paste("beta",1:d,sep="_"))
	
	if( length(h) == 0 )
	{
		h <- bw.SJ(model$resid)
	} 
	
	BBgetweights.out <- BBgetweights(n,l,B)	
	MBB.weights <- BBgetweights.out$MBBweights
	ETBB.weights <- BBgetweights.out$ETBBweights
	m.l <- BBgetweights.out$m_l

	pi.tilde.MBB <- get.pi.tilde.MBB(n,l)
	pi.tilde.ETBB <- get.pi.tilde.ETBB(n,l)

	# get appropriate centering vectors for different bootstrap methods:
	
	MBB.beta.tilde <- rq(	Y ~ X[,-1],
							tau = tau,
							weights = pi.tilde.MBB)$coef

	ETBB.beta.tilde <- rq(	Y ~ X[,-1],
							tau = tau,
							weights = pi.tilde.ETBB)$coef
								
	SMBB.beta.tilde <- nlm(	f = S.tilde,
							p = beta.hat,
							Y = Y,
							X = X,
							pi.tilde = pi.tilde.MBB,
							h = h,
							tau = tau,
							iterlim = 200)$estimate
										
	SETBB.beta.tilde <- nlm(f = S.tilde,
							p = beta.hat,
							Y = Y,
							X = X,
							pi.tilde = pi.tilde.ETBB,
							h = h,
							tau = tau,
							iterlim = 200)$estimate

	MBB.beta.hat <- ETBB.beta.hat <- SMBB.beta.hat <- SETBB.beta.hat <- matrix(NA,B,d+1)
	
	for(b in 1:B)
	{	
		
		MBB.beta.hat[b,] <- rq(	Y ~ X[,-1],
								tau = tau,
								weights = MBB.weights[,b])$coef
		
		ETBB.beta.hat[b,] <- rq(Y ~ X[,-1],
								tau = tau,
								weights = ETBB.weights[,b])$coef
			
		Y.smooth <- Y + rnorm(n,0,h)
		X.smooth <- X + cbind(rep(0,n),matrix(rnorm(n*d,0,h),n,d))
		
		SMBB.beta.hat[b,] <- rq(Y.smooth ~ X.smooth[,-1],
								tau = tau,
								weights = MBB.weights[,b])$coef

		SETBB.beta.hat[b,] <- rq(Y.smooth ~ X.smooth[,-1],
								tau = tau,
								weights = ETBB.weights[,b])$coef			

	}

	MBB.cov.est <- cov(MBB.beta.hat)
	ETBB.cov.est <- m.l * cov(ETBB.beta.hat)

	SMBB.cov.est <- cov(SMBB.beta.hat)
	SETBB.cov.est <- m.l * cov(SETBB.beta.hat)
	
	MBB.pivot <- sqrt(n) * ( MBB.beta.hat - matrix(MBB.beta.tilde,B,d+1,byrow=TRUE))
	ETBB.pivot <- sqrt(m.l) * sqrt(n) * ( ETBB.beta.hat - matrix(ETBB.beta.tilde,B,d+1,byrow=TRUE))
	SMBB.pivot <- sqrt(n) * ( SMBB.beta.hat - matrix(SMBB.beta.tilde,B,d+1,byrow=TRUE))
	SETBB.pivot <- sqrt(m.l) * sqrt(n) * ( SETBB.beta.hat - matrix(SETBB.beta.tilde,B,d+1,byrow=TRUE))	
	
	MBB.loci <- beta.hat - apply(MBB.pivot,2,quantile,1-alpha/2)/sqrt(n)
	MBB.upci <- beta.hat - apply(MBB.pivot,2,quantile,alpha/2)/sqrt(n)
	MBB.confint <- cbind(MBB.loci,MBB.upci)
	
	ETBB.loci <- beta.hat - apply(ETBB.pivot,2,quantile,1-alpha/2)/sqrt(n)
	ETBB.upci <- beta.hat - apply(ETBB.pivot,2,quantile,alpha/2)/sqrt(n)
	ETBB.confint <- cbind(ETBB.loci,ETBB.upci)
	
	SMBB.loci <- beta.hat - apply(SMBB.pivot,2,quantile,1-alpha/2)/sqrt(n)
	SMBB.upci <- beta.hat - apply(SMBB.pivot,2,quantile,alpha/2)/sqrt(n)
	SMBB.confint <- cbind(SMBB.loci,SMBB.upci)
	
	SETBB.loci <- beta.hat - apply(SETBB.pivot,2,quantile,1-alpha/2)/sqrt(n)
	SETBB.upci <- beta.hat - apply(SETBB.pivot,2,quantile,alpha/2)/sqrt(n)
	SETBB.confint <- cbind(SETBB.loci,SETBB.upci)
	
	output <- list( beta.hat = beta.hat,
	                MBB.pivot = MBB.pivot,
					        ETBB.pivot = MBB.pivot,
        					SMBB.pivot = SMBB.pivot,
        					SETBB.pivot = SETBB.pivot,
        					MBB.confint = MBB.confint,
        					ETBB.confint = ETBB.confint,
        					SMBB.confint = SMBB.confint,
        					SETBB.confint = SETBB.confint,
        					MBB.cov.est = MBB.cov.est,
        					ETBB.cov.est = ETBB.cov.est,
        					SMBB.cov.est = SMBB.cov.est,
        					SETBB.cov.est = SETBB.cov.est
					      )
	
	output$call <- match.call()
	
	class(output) <- "QregBB"
	
	return(output)
	
}


#' Print method for class \code{QregBB}
#' @export
#' @noRd
print.QregBB <- function(x,...){
  
  cat("\nCall:\n",
      paste(deparse(x$call), sep="\n", collapse = "\n"), "\n", sep = "")
  
  cat("\nCoefficients:\n")
  
  table <- as.table(cbind(round(x$beta.hat,5),
                          round(sqrt(diag(x$MBB.cov.est)),5),
                          round(sqrt(diag(x$SMBB.cov.est)),5),
                          round(sqrt(diag(x$ETBB.cov.est)),5),
                          round(sqrt(diag(x$SETBB.cov.est)),5)))
  
  colnames(table) <- c("Estimate","SE (MBB)","SE (SMBB)","SE (ETBB)","SE (SETBB)")
  rownames(table) <- names(x$beta.hat)
  print(table)
  
  cat("\nConfidence intervals:\n")
  
  table <- as.table(cbind(round(x$beta.hat,5),
                          round(x$MBB.confint,5),
                          round(x$SMBB.confint,5),
                          round(x$ETBB.confint,5),
                          round(x$SETBB.confint,5)))
  
  colnames(table) <- c("Estimate",
                       "lower (MBB)","upper (MBB)",
                       "lower (SMBB)","upper (SMBB)",
                       "lower (ETBB)","upper (ETBB)",
                       "lower (SETBB)","upper (SETBB)")
  
  rownames(table) <- names(x$beta.hat)
  
  print(table)
  
}

#' A function asking how many members of I are members of J
#' 
#' @param I a vector
#' @param J a vector
#' @return the number of entries in \code{I} which are also in \code{J}.
#' @noRd
IinJ <- function(I,J){sum(I %in% J) == 0}

#' Computes D.star from paper
#'
#' @param Y the vector of response values.
#' @param X the design matrix (including a column of ones for the intercept).
#' @param beta parameter vector at which to compute the objective function.
#' @param pi.star weights retrieved from the \code{BBgetweights} function.
#' @param tau the quantile of interest.
#' @returns the value of \eqn{\tilde D_n^*} from the paper.
#' @noRd
D.n.star <- function(Y,X,beta,tau,pi.star)
{
	
	apply( - pi.star * X * ( tau - as.numeric( Y - X %*% beta <= 0)),2,sum)
	
}

#' Chooses block sizes for MBB, ETBB, SMBB, and SETBB via the NPPI for quantile regression
#'
#' @param Y the vector of response values.
#' @param X the design matrix (including a column of ones for the intercept).
#' @param tau the quantile of interest.
#' @param min.in.JAB the minimum number of Monte-Carlos draws desired in each jackknife draw
#' @return Returns a list of the NPPI-selected block sizes for the MBB, SMBB, ETBB, and SETBB.
#' 
#' @details This function is based on the nonparametric plug-in (NPPI) method discussed in Lahiri (2003), which makes use of the jackknife-after-bootstrap (JAB).
#'
#' @references 
#' 
#' Gregory, K. B., Lahiri, S. N., & Nordman, D. J. (2018). A smooth block bootstrap for quantile regression with time series. \emph{The Annals of Statistics}, 46(3), 1138-1166.
#' 
#' Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer, New York.
#'
#' @examples
#' # generate some data and use NPPI to choose block sizes for MBB, SMBB, ETBB, and SETBB.
#' n <- 50
#' X1 <- arima.sim(model=list(ar=c(.7,.1)),n)
#' X2 <- arima.sim(model=list(ar=c(.2,.1)),n)
#' e <- arima.sim(model=list(ar=c(.7,.1)),n)
#' Y <- X1 + e
#' X <- cbind(rep(1,n),X1,X2)
#' 
#' blksize.out <- getNPPIblksizesQR(Y,X,tau=.5)
#' blksize.out
#' @export
getNPPIblksizesQR <- function(Y,X,tau,min.in.JAB=100)
{
		
	n <- length(Y)
	d <- ncol(X) - 1
	
	l.1 <- round( n^(1/5)) # block size to use for estimation of B.hat and v.hat: (7.63) pg. 196 of Lahiri (2003) 
	m <- floor( n^(1/3)*l.1^(2/3)) # block length for JAB--number of contiguous blocks to remove: (7.63) pg. 196 of Lahiri (2003)
	l.2 <- 2*l.1

	B <- round(min.in.JAB*( 1 - m/(n - l.1 + 1))^(-floor(n/l.1))) 

	pi.tilde.MBB.l.1 <- get.pi.tilde.MBB(n,l.1)	
	pi.tilde.MBB.l.2 <- get.pi.tilde.MBB(n,l.2)	
	pi.tilde.ETBB.l.1 <- get.pi.tilde.ETBB(n,l.1)
	pi.tilde.ETBB.l.2 <- get.pi.tilde.ETBB(n,l.2)
	
	model <- rq(Y ~ X+0, tau = tau)
	beta.hat <- model$coef
	
	sj.bw <- bw.SJ(model$resid)

	perturbations <- array(rnorm(n = n * (d+1) * B),dim=c(n,d+1,B))

	MBB.boot.D.n.l.1 <- ETBB.boot.D.n.l.1 <- matrix(NA,B,d+1)	
	SMBB.boot.D.n.l.1  <- SETBB.boot.D.n.l.1<- matrix(NA,B,d+1)

	MBB.boot.D.n.l.2 <- ETBB.boot.D.n.l.2 <- matrix(NA,B,d+1)
	SMBB.boot.D.n.l.2  <- SETBB.boot.D.n.l.2<- matrix(NA,B,d+1)

	# l.1
				
	BBweights.out.l.1 <- BBgetweights( n, l.1 ,B= B)
	pi.star.MBB.l.1 <- BBweights.out.l.1$MBBweights
	pi.star.ETBB.l.1 <- BBweights.out.l.1$ETBBweights
	m.l.l.1 <- BBweights.out.l.1$m_l
	
	beta.tilde.MBB.l.1 <- rq(	Y ~ X+0,
								tau =  tau,
								weights = pi.tilde.MBB.l.1)$coef
	
	beta.tilde.ETBB.l.1 <- rq(	Y ~ X+0,
								tau =  tau,
								weights = pi.tilde.ETBB.l.1)$coef
								
	beta.tilde.SMBB.l.1 <- nlm(	f = S.tilde,
								p = beta.hat,
								Y = Y,
								X = X,
								pi.tilde = pi.tilde.MBB.l.1,
								h = sj.bw,
								tau =  tau,
								iterlim = 200)$estimate
											
	beta.tilde.SETBB.l.1 <- nlm(	f = S.tilde,
								p = beta.hat,
								Y = Y,
								X = X,
								pi.tilde = pi.tilde.ETBB.l.1,
								h = sj.bw,
								tau =  tau,
								iterlim = 200)$estimate
								
	# l.2
		
	BBweights.out.l.2 <- BBgetweights( n, l.2 ,B= B)
	pi.star.MBB.l.2 <- BBweights.out.l.2$MBBweights
	pi.star.ETBB.l.2 <- BBweights.out.l.2$ETBBweights
	m.l.l.2 <- BBweights.out.l.2$m_l							
								
	beta.tilde.MBB.l.2 <- rq(	Y ~ X+0,
								tau =  tau,
								weights = pi.tilde.MBB.l.2)$coef
	
	beta.tilde.ETBB.l.2 <- rq(	Y ~ X+0,
								tau =  tau,
								weights = pi.tilde.ETBB.l.2)$coef
								
	beta.tilde.SMBB.l.2 <- nlm(	f = S.tilde,
								p = beta.hat,
								Y = Y,
								X = X,
								pi.tilde = pi.tilde.MBB.l.2,
								h = sj.bw,
								tau =  tau,
								iterlim = 200)$estimate
											
	beta.tilde.SETBB.l.2 <- nlm(	f = S.tilde,
								p = beta.hat,
								Y = Y,
								X = X,
								pi.tilde = pi.tilde.ETBB.l.2,
								h = sj.bw,
								tau =  tau,
								iterlim = 200)$estimate
	
		for(b in 1:B)
		{
						
			Y.smooth <- Y + rnorm(n,0,sj.bw)
			X.smooth <- X + cbind(rep(0,n),matrix(rnorm(n*d,0,sj.bw),n,d))
		
			# l.1
																								
			MBB.boot.D.n.l.1[b,] <- D.n.star(Y,X,beta.tilde.MBB.l.1,tau= tau,pi.star.MBB.l.1[b])		
			
			ETBB.boot.D.n.l.1[b,] <- D.n.star(Y,X,beta.tilde.ETBB.l.1,tau= tau,pi.star.ETBB.l.1[b])						
																					
			SMBB.boot.D.n.l.1[b,] <- D.n.star(Y,X,beta.tilde.SMBB.l.1,tau= tau,pi.star.MBB.l.1[b])		
			
			SETBB.boot.D.n.l.1[b,] <- D.n.star(Y,X,beta.tilde.SETBB.l.1,tau= tau,pi.star.ETBB.l.1[b])		
			
			# l.2 
												
			MBB.boot.D.n.l.2[b,] <- D.n.star(Y,X,beta.tilde.MBB.l.2,tau= tau,pi.star.MBB.l.2[b])		
			
			ETBB.boot.D.n.l.2[b,] <- D.n.star(Y,X,beta.tilde.ETBB.l.2,tau= tau,pi.star.ETBB.l.2[b])						
																					
			SMBB.boot.D.n.l.2[b,] <- D.n.star(Y,X,beta.tilde.SMBB.l.2,tau= tau,pi.star.MBB.l.2[b])		
			
			SETBB.boot.D.n.l.2[b,] <- D.n.star(Y,X,beta.tilde.SETBB.l.2,tau= tau,pi.star.ETBB.l.2[b])		
															
		}
		


# Choose block sizes with NPPI - okay, we will pretend that we are trying to estimate the trace of Sigma.
# There is a function of X with asymptotic variance equal to this trace...


			
	M <- ( n-l.1+1) - m + 1 # sample size after removing m blocks of length l.1
	N <-  n - l.1 + 1
	
	MBB.JAB.phi.hat <- ETBB.JAB.phi.hat <- SMBB.JAB.phi.hat <- SETBB.JAB.phi.hat <- numeric()
	
	for(k in 1:M)
	{

		which.for.k <- apply(BBweights.out.l.1$blkpoints, 2, IinJ , I = c(k:(k+m-1)) )  
		MBB.JAB.phi.hat[k] <-  sum(diag(cov(MBB.boot.D.n.l.1[which.for.k,]))) 	# estimated variance of sum(D.n) from MBB bootstraps which.for.k, block size l.1
		ETBB.JAB.phi.hat[k] <-  sum(diag(cov(ETBB.boot.D.n.l.1[which.for.k,])))	# estimated variance of sum(D.n) from ETBB bootstraps which.for.k, block size l.1
		SMBB.JAB.phi.hat[k] <-  sum(diag(cov(SMBB.boot.D.n.l.1[which.for.k,]))) # estimated variance of sum(D.n) from MBB bootstraps which.for.k, block size l.1
		SETBB.JAB.phi.hat[k] <-  sum(diag(cov(SETBB.boot.D.n.l.1[which.for.k,])))# estimated variance of sum(D.n) from ETBB bootstraps which.for.k, block size l.1			
	}
	
	MBB.phi <- sum(diag(cov(MBB.boot.D.n.l.1))) # estimated trace of cov(D.n) from all MBB bootstraps, block size l.1
	MBB.JAB.phi.tilde <- (1/m)*( N * MBB.phi - (N - m) * MBB.JAB.phi.hat) # pseudo values, cf. pg. 192 Lahiri (2003)
	var.hat.JAB.MBB <- m/(N - m) * mean((MBB.JAB.phi.tilde - MBB.phi )^2) # formula (7.53), pg. 192 of Lahiri (2003)
	v.hat.JAB.MBB <-   n/l.1 * var.hat.JAB.MBB # pg. 195 of Lahiri (2003)
	B.hat.MBB <- 2*l.1*(sum(diag(cov(MBB.boot.D.n.l.1))) - sum(diag(cov(MBB.boot.D.n.l.2)))) # pg. 195 of Lahiri (2003)
	l.opt.MBB <- min(round(n/2),max(1,round(  n^(1/3) * ( 2 * B.hat.MBB^2 / v.hat.JAB.MBB )^(1/3) ))) # eq. (7.57), pg. 194 of Lahiri(2003)
		
	ETBB.phi <- sum(diag(cov(ETBB.boot.D.n.l.1))) # estimated trace of cov(D.n) from all ETBB bootstraps, block size l.1
	ETBB.JAB.phi.tilde <- (1/m)*( N * ETBB.phi - (N - m) * ETBB.JAB.phi.hat) # pseudo values, cf. pg. 192 Lahiri (2003)
	var.hat.JAB.ETBB <- m/(N - m) * mean((ETBB.JAB.phi.tilde - ETBB.phi )^2) # formula (7.53), pg. 192 of Lahiri (2003)
	v.hat.JAB.ETBB <-   n/l.1 * var.hat.JAB.ETBB # pg. 195 of Lahiri (2003)
	B.hat.ETBB <- (4/3)*l.1^2*(sum(diag(cov(ETBB.boot.D.n.l.1))) - sum(diag(cov(ETBB.boot.D.n.l.2)))) # pg. 195 of Lahiri (2003)
	l.opt.ETBB <- min(round(n/2),max(1,round(  n^(1/5) * ( 4 * B.hat.ETBB^2 / v.hat.JAB.ETBB )^(1/5)) )) # different bias order for tapered blocks
			
	SMBB.phi <- sum(diag(cov(SMBB.boot.D.n.l.1))) # estimated trace of cov(D.n) from all SMBB bootstraps, block size l.1
	SMBB.JAB.phi.tilde <- (1/m)*( N * SMBB.phi - (N - m) * SMBB.JAB.phi.hat) # pseudo values, cf. pg. 192 Lahiri (2003)
	var.hat.JAB.SMBB <- m/(N - m) * mean((SMBB.JAB.phi.tilde - SMBB.phi )^2) # formula (7.53), pg. 192 of Lahiri (2003)
	v.hat.JAB.SMBB <-   n/l.1 * var.hat.JAB.SMBB # pg. 195 of Lahiri (2003)
	B.hat.SMBB <- 2*l.1*(sum(diag(cov(SMBB.boot.D.n.l.1))) - sum(diag(cov(SMBB.boot.D.n.l.2)))) # pg. 195 of Lahiri (2003)
	l.opt.SMBB <- min(round(n/2),max(1,round(  n^(1/3) * ( 2 * B.hat.SMBB^2 / v.hat.JAB.SMBB )^(1/3) ))) # eq. (7.57), pg. 194 of Lahiri(2003)

	SETBB.phi <- sum(diag(cov(SETBB.boot.D.n.l.1))) # estimated trace of cov(D.n) from all SETBB bootstraps, block size l.1
	SETBB.JAB.phi.tilde <- (1/m)*( N * SETBB.phi - (N - m) * SETBB.JAB.phi.hat) # pseudo values, cf. pg. 192 Lahiri (2003)
	var.hat.JAB.SETBB <- m/(N - m) * mean((SETBB.JAB.phi.tilde - SETBB.phi )^2) # formula (7.53), pg. 192 of Lahiri (2003)
	v.hat.JAB.SETBB <-   n/l.1 * var.hat.JAB.SETBB # pg. 195 of Lahiri (2003)
	B.hat.SETBB <- (4/3)*l.1^2*(sum(diag(cov(SETBB.boot.D.n.l.1))) - sum(diag(cov(SETBB.boot.D.n.l.2)))) # pg. 195 of Lahiri (2003)
	l.opt.SETBB <- min(round(n/2),max(1,round(  n^(1/5) * ( 4 * B.hat.SETBB^2 / v.hat.JAB.SETBB )^(1/5)) )) # different bias order for tapered blocks

out <- list(l.opt.MBB  = l.opt.MBB ,
			l.opt.ETBB = l.opt.ETBB,
			l.opt.SMBB = l.opt.SMBB,
			l.opt.SETBB = l.opt.SETBB)
			
return(out)
					
}

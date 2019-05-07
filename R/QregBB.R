#' Retrieve weights from from Monte-Carlo block samples for the ETBB
#'
#' @param n length of observed time series.
#' @param l block length.
#' @param B the number of Monte-Carlo draws.
#' @param seed a random number generating seed (optional).
#' @param c a constant controlling the tapering of weight assignments over the blocks.
#' @return a list containing an \code{n} by \code{B} matrix containing the the mass placed on the \code{n} observations by the ETBB empirical distribution in each of the \code{B} Monte-Carlo draws. Another object in the list gives more info.
ETBBgetweights<-function(n,l,B,seed=NA,c=.43){

	if(is.na(seed)) seed<-floor(runif(1)*100000)

	ETBBweights<-.C("ETBBgetweights",
	n = as.integer(n),
	l = as.integer(l),
	B = as.integer(B),
	blkpts = integer(floor(n/l+1)),
	f_tilde = double(B*n+l),
	as.integer(seed),
	weights = double(l),
	as.double(c))		

	weights <- matrix( ETBBweights$f_tilde[1:(n*B)],n,B)

	output <-  list(ETBBweights = weights,
			other = ETBBweights)

	return(output)
}


#' Retrieve weights from from Monte-Carlo block samples for the MBB
#'
#' @param n length of observed time series.
#' @param l block length.
#' @param B the number of Monte-Carlo draws.
#' @param seed a random number generating seed (optional).
#' @return a list containing an \code{n} by \code{B} matrix containing the the mass placed on the \code{n} observations by the MBB empirical distribution in each of the \code{B} Monte-Carlo draws. Another object in the list gives more info.
MBBgetweights<-function(n,l,B,seed=NA){

	if(is.na(seed)) seed<-floor(runif(1)*100000)

	MBBweights<-.C("MBBgetweights",
	n = as.integer(n),
	l = as.integer(l),
	B = as.integer(B),
	blkpts = integer(floor(n/l+1)),
	f_tilde = double(B*n+l),
	as.integer(seed))		

	weights <- matrix( MBBweights$f_tilde[1:(n*B)],n,B)

	output <-  list(MBBweights = weights,
					other = MBBweights)

	return(output)
}

#' Retrieve weights from from Monte-Carlo block samples for the MBB and ETBB
#'
#' @param n length of observed time series.
#' @param l block length.
#' @param B the number of Monte-Carlo draws.
#' @param c a constant controlling the tapering of weight assignments over the blocks.
#' @param seed a random number generating seed (optional).
#' @return a list containing \code{n} by \code{B} matrices containing the the mass placed on the \code{n} observations by the MBB and ETBB empirical distributions in each of the \code{B} Monte-Carlo draws, the sampled indices for the block starting points, and the scalar \eqn{m_l} corresponding to the weight function and the block length.
AllBBgetweights <- function(n,l,B,seed=NA,c=.43){

	if(is.na(seed)) seed <- floor(runif(1)*100000)

	AllBBweights<-.C("AllBBgetweights",
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

	ETBBweights <- matrix( AllBBweights$f_ETBB[1:(n*B)],n,B)
	MBBweights <- matrix( AllBBweights$f_MBB[1:(n*B)],n,B)
	blkpoints <-  matrix( AllBBweights$blkpts,floor(n/l),B) + 1

	output <-  list(ETBBweights = ETBBweights,
					MBBweights = MBBweights,
					blkpoints = blkpoints,
					m_l = AllBBweights$m_l )


	return(output)
}


#' Function defining the tapered assignment of weights across the ETBB blocks
#'
#' @param u a real number.
#' @param c a constant controlling the tapering.
#' @return the value of the tapering function.
omega <- function( u,c=.43){ ( u <= c) * u/c + ((u > c) & (u < 1 - c)) + (u >= 1-c) *(1 - u)/c}

#' Computes expected value of the bootstrap weights assigned to each time point by the ETBB
#'
#' @param n the length of the observed time series.
#' @param l the block length.
#' @param c a number controlling the tapering.
#' @return a vector of length \code{n} containing the expected values of the weights assigned to each time point by the ETBB.
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

S.tilde <- function(beta,Y,X,pi.tilde,h,tau)
{
	
	var.eta <- h^2*(1 + sum(beta[-1]^2))
	e.hat <- Y - X %*% beta
	Phi.e.hat <- pnorm(e.hat,0,sqrt(var.eta))
	phi.e.hat <- dnorm(e.hat,0,sd=sqrt(var.eta))
	S.tilde <- sum(pi.tilde*( e.hat*(tau-(1-Phi.e.hat)) + var.eta*phi.e.hat) )
	
	return(S.tilde)

}

#' Computes bootstrap expected value of the quantile regression objective function
#'
#' @param Y the vector of response values.
#' @param X the design matrix (including a column of ones for the intercept).
#' @param tau the quantile of interest.
#' @param l block size.
#' @param B the number of Monte Carlo bootstrap samples to draw.
#' @param h a scalar bandwidth (bandwidth matrix is \code{h} times identity).
#' @param alpha a significance level to which the returned confidence intervals will correspond.
#' @return a list containing for the MBB, SMBB, ETBB, and SETBB the set of Monte Carlo draws of the pivot quantity \eqn{\sqrt{n}(\hat \beta^*_n - \tilde \beta_n)}, confidence intervals for each component of \eqn{\beta} corresponding to the specified confidence level, and estimates of the asymptotic covariance matrix of the pivot quantity \eqn{\sqrt{n}(\hat \beta_n - \beta)}.
QregBB <- function(Y,X,tau,l,B=500,h=NULL,alpha=0.05)
{
	
	n <- length(Y)
	d <- ncol(X) - 1

	model <- rq(Y ~ X[,-1], tau = tau)
	beta.hat <- model$coef
	if( length(h)==0)
	{
		h <- bw.SJ(model$resid)
	} 
	
	AllBBgetweights.out <- AllBBgetweights(n,l,B)	
	MBB.weights <- AllBBgetweights.out$MBBweights
	ETBB.weights <- AllBBgetweights.out$ETBBweights
	m.l <- AllBBgetweights.out$m_l

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
	
	output <- list( MBB.pivot = MBB.pivot,
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
					
}
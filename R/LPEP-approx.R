
#### Load libraries #####

library(BayesLogit)
library(mvtnorm)
library(robustbase)
library(WoodburyMatrix)
library(testit)
library(detectseparation)
library(ncvreg)



#### Helper functions ####

#### Draws from reflective normal distribution with left reflection ####
# boundary a and given mean and sd

#' Random generation for the reflective Gaussian distribution
#'
#' @param n number of observations. If length(n)>1, the length is taken to be the number required
#' @param a left reflection boundary. Default is 0
#' @param mean mean
#' @param sd standard deviation
#'
#' @return rrefnorm generates random deviates from the specified reflective gausisan distribution
#' @export
#'
#' @examples
rrefnorm <- function(n,a=0, mean=0, sd=1){
  eps <- rnorm(n)
  g=a+abs(mean+sd*eps-a)
  return(g)
}

#### Calculates density of reflective normal distribution  with left reflection boundary ####
#  a and given mean and sd
#' Density calculation for the reflective Gaussian distribution
#'
#' @param x vector of quantiles
#' @param a left reflection boundary default is 0
#' @param mean mean
#' @param sd standard deviation
#' @param logarithm logical; if TRUE, density evaluated on log scale
#'
#' @return returns the density of the specified reflective Gaussian distribution at x.
#' @export
#'
#' @examples
drefnorm <- function(x,a=0,mean=0,sd=1,logarithm=FALSE){
  if(x>=a){
    val <- dnorm(x,mean=mean,sd=sd)+dnorm(x,mean=2*a-mean, sd=sd)
  }else{
    val <- 0
  }
  if(logarithm==TRUE){return(log(val))}
  return(val)
}



#### Calculates density of hyper-g/n prior with specified alpha and n ####
dhypergn <- function(g, alpha=4,n,logarithm=FALSE){
  val <- (alpha-2)/(2*n)*(1/(1+g/n))^(alpha/2)
  if(logarithm==TRUE){
    val <- log(alpha-2)-log(2*n)-alpha/2*log(1+g/n)
  }
  return (val)
}

#### Draws from hyper-g/n distribution with specified alpha and n ####
rhypergn <- function(nsamp, alpha=4,n){
  u <- rbeta(nsamp, alpha/2-1,1)
  return (n*u/(1-u))
}


#### Calculates density of a robust prior with specified a,b, rho and n ####
drobust <- function(g, a=0.5,b=1,n,rho,logarithm=FALSE){
  if(g>rho*(b+n)-b){
    val <- a*(rho*(b+n))^a*(g+b)^(-a-1)
  }else{ val <- 0}

  if(logarithm==TRUE){
    if(g>rho*(b+n)-b){
      val <- log(a)+a*log(rho*(b+n))-(a+1)*log(g+b)
    }else{
      val <- -Inf
    }
  }
  return (val)
}

#### Function to evaluate MH acceptance probability of the proposed delta ####

# choices of delta allowed are gamma, hyper, hyper-g/n and robust #
delta.acceptance.prob <- function(delta.new,delta.old,bmat, y.star,xmat,gam,hyper.type,hyper.param,exact.mixture.g,seb.held){
  n <- length(y.star)
  if(sum(gam)==0){
    mod1 <- glm(y.star~1,family = binomial(link = "logit"))
  }else{
    mod1 <- glm(y.star~xmat[,as.logical(gam)],family = binomial(link = "logit"))
  }


  if(hyper.type=="gamma") {prior.d <-
    dgamma(delta.new/n, shape = hyper.param, rate = hyper.param, log = TRUE) -
    dgamma(delta.old/n, shape = hyper.param, rate = hyper.param, log = TRUE)}
  if (hyper.type=="hyper-g/n") {prior.d <-
    dhypergn(delta.new,alpha = hyper.param,n=n,logarithm = TRUE)-
    dhypergn(delta.old,alpha = hyper.param,n=n,logarithm = TRUE)}
  if (hyper.type=="hyper-g") {prior.d <-
    dhypergn(delta.new,alpha = hyper.param,n=1,logarithm = TRUE)-
    dhypergn(delta.old,alpha = hyper.param,n=1,logarithm = TRUE)}
  if(hyper.type=="robust"){prior.d <-
    drobust(delta.new, n=n, rho=1/(sum(gam)+1),logarithm = TRUE)-
    drobust(delta.old, n=n, rho=1/(sum(gam)+1),logarithm = TRUE)}
  if(exact.mixture.g==TRUE){
    beta.star <- rep(0, length(mod1$coefficients))
    inf.mat <- vcov(mod1)
  }else if(seb.held==TRUE){
    beta.star <- rep(0, length(mod1$coefficients))
    Xgam <- model.matrix(mod1)
    inf.mat <- 4*solve(t(Xgam)%*% Xgam)
  }else{
    beta.star <- mod1$coefficients
    inf.mat <- vcov(mod1)
  }


  log.num <- prior.d + dmvnorm(t(bmat), mean = beta.star, sigma = delta.new*inf.mat,log=TRUE)
  log.den <- dmvnorm(t(bmat), mean = beta.star, sigma = delta.old*inf.mat,log=TRUE)

  a.p.delta <- min(exp(log.num-log.den),1)

  return (a.p.delta)
}

#### Proposal function for imaginary samples ####
proposal.ystar <- function(y.star, pi.star, n_i, d=5){
  n <- length(y.star)
  prob1 <- c(0.5,0.2,0.15,0.10,0.05)
  prob2 <- c(0.7,0.3)

  s <- resample(1:2,1, prob = prob2)
  # Use technique of George mchlloch stat sinaca paper  eq 46 if s==1
  # Local move
  if(s==1){
    d0 <- resample(1:d,1, prob=prob1)
    ind <- resample(1:n,d0)
    y.star.cand <- y.star
    y.star.cand[ind]=as.numeric(!y.star.cand[ind])
    # Global move
  }else if(s==2){


    y.star.cand <- rbinom(n, n_i,  pi.star)

  }
  return(list(y.star.cand=y.star.cand,s=s))
}

#### Function to evaluate MH acceptance probability of the proposed imaginary sample ####
ystar.acceptance.prob <- function(y.new,y.old,xmat,gam, bmat,delta,n_i,pi.star,local,mystar.mod){
  n <- length(y.new)
  new.sum <- sum(y.new)
  old.sum <- sum(y.old)
  if(sum(gam)==0){
    mod.new <- glm(y.new~1,family = binomial(link = "logit"))
    mod.old <- glm(y.old~1,family = binomial(link = "logit"))
  }else{
    mod.new <- glm(y.new~xmat[, as.logical(gam)],family = binomial(link = "logit"))
    mod.old <- glm(y.old~xmat[, as.logical(gam)],family = binomial(link = "logit"))
  }

  if(is.null(mystar.mod)){
    dmystar.new <- lgamma(new.sum+0.5)+lgamma(n-new.sum+0.5)
    dmystar.old <- lgamma(old.sum+0.5)+lgamma(n-old.sum+0.5)
  }else{
    prob.ystar <- predict(mystar.mod,xmat,type="response")
    dmystar.new <- sum(dbinom(y.new,1,prob.ystar,log = TRUE))
    dmystar.old <-sum(dbinom(y.old,1,prob.ystar,log=TRUE))
  }

  new.log.bgam.dens <- dmvnorm(t(bmat), mean = mod.new$coefficients, sigma = delta*vcov(mod.new),log=TRUE)
  old.log.bgam.dens <- dmvnorm(t(bmat), mean = mod.old$coefficients, sigma = delta*vcov(mod.old),log=TRUE)

  log.num <- dmystar.new + new.log.bgam.dens
  +ifelse(local==2,dbinom(y.old,n_i,prob=pi.star,log=TRUE),0)
  log.den <- dmystar.old + old.log.bgam.dens
  +ifelse(local==2,dbinom(y.new,n_i,prob=pi.star,log=TRUE),0)

  a.p <- min(exp(log.num-log.den),1)

  return (a.p)

}



resample <- function(x, ...) x[sample.int(length(x), ...)]


#### Proposal function for gamma variable ####
proposal.gamma <- function(gam, d=4){
  p <- length(gam)
  prob1 <- c(0.6,0.2,0.15,0.05)
  prob2 <- c(0.9,0.1)

  s <- resample(1:2,1, prob = prob2)
  # Use technique of George mchlloch stat sinaca paper  eq 46 if s==1

  if(s==1){
    if(p<d){
      d=p
      prob1 <- prob1[1:d]/sum(prob1[1:d])
    }
    d0 <- resample(1:d,1, prob=prob1)
    ind <- resample(1:p,d0)
    gam[ind]=!gam[ind]


    # swap one entry from active set with one entry from non-active set randomly
  }else if(s==2){
    if(all(gam==1)){
      ind <- resample(1:p,1)
      gam[ind]=0
    }else if(all(gam==0)){
      ind <- resample(1:p,1)
      gam[ind]=1
    }else{
      gam1 <- which(gam==1)
      gam0 <- which(gam==0)
      ind1 <- resample(gam1,1)
      ind0 <- resample(gam0,1)
      gam[ind1] <- 0
      gam[ind0] <- 1
    }
  }
  return(gam)
}


#### Expit function ####
expit <- function(x){
  return(1/(1+exp(-x)))
}

#### Bayesian Inference for logistic models with Laplace PEP prior methodology ####

#' Bayesian Inference for logistic models with Laplace PEP prior methodology

#'
#' @param x Matrix of covariates, dimension is n*p; Intercept does not need to be included
#' @param y Response, a n*1 binary vector
#' @param burn Number of burn-in MCMC iterations. Default is 1000
#' @param nmc Number of MCMC iterations to be saved post burn-in. Default is 5000
#' @param model.prior Family of prior distribution on the models. Currently, two choices allowed: Uniform and beta-binomial; Default is beta-binomial.
#' @param hyper Logical; True if hyper prior on delta is specified. Default is TRUE. If FALSE, UIP version is fit.
#' @param hyper.type If hyper= TRUE, prior type on delta can be specified here.
#' Currently, four choices are offered: "gamma", "hyper-g", "hyper-g/n" and "robust". Under
#' gamma prior, delta/n follows a gamma with shape and rate parameter equal to hyper.param.
#' Under hyper-g and hyper-g/n parameter, value of a parameter is specified by hyper.param.
#' For robust, recommended choices in Bayarri et Al 2012 for
#' a=0.5, b=1 and rho=1/(pgam +1) are used.
#'
#' @param hyper.param Value of hyper-parameter for prior on delta if any.
#' If hyper.type="gamma", shape=rate= hyper.param value and hence hyper.param>0.
#' If hyper.type="hyper-g", or hyper.type="hyper-g/n", value of a=hyper.param and hence hyper.param >2
#' with recommended values to be 3 or 4.
#' @param exact.mixture.g Logical; If True, exact version of Li and Clyde (2018) mixture of g-priors
#' is implemented using Polya-gamma augmentation instead of Laplace approximation. Default is False.
#' @param seb.held Logical; If True, Bove et Al (2011) prior is implemented
#'
#' @return Laplace.pep returns a fit object. This object contains the following components:
#'
#' @export
#'
#' @examples
Laplace.pep.approx <- function(x,y,
                        burn=1000,
                        nmc=5000,
                        model.prior="beta-binomial",
                        hyper="TRUE",
                        hyper.type="hyper-g/n",
                        hyper.param=NULL,
                        mystar=NULL,
                        exact.mixture.g=FALSE,
                        seb.held=FALSE,
                        sample.beta=TRUE){

  n <- length(y)
  p <- ncol(x)
  # sd for delta proposal i.e. reflective gaussian sd
  tau=n/2

  # create matrix to save variables
  GammaSave = matrix(NA, nmc, p)
  BetaSave = matrix(NA, nmc, p+1)
  OmegaSave = matrix(NA, nmc, n)
  ystarSave = matrix(NA, nmc, n)
  deltaSave = matrix(NA,nmc, 1)
  timemat <- matrix(NA, nmc+burn, 5)
  #sep.time <- matrix(NA, nmc+burn,1)
  # Intialize parameters
  gam = rep(0,p)
  ful.mod <- glm(y~1,family = binomial(link = "logit"))
  b = c(ful.mod$coefficients,rep(0,p))
  n_i = rep(1,n)
  # omega = unlist(lapply(n_i, rpg,num=1,z=0)) # prior is PG(n_i, 0) where n_i=1 for binary logit

  if(exact.mixture.g==TRUE){
    ystar=y
  }else{
    ystar = rbinom(n,1,0.5)
  }
  delta = n

  if(is.null(mystar)){
    mystar.mod <- NULL
  }else if(mystar=="SCAD"){
    mystar.mod <-  cv.ncvreg(x,y,family = "binomial", penalty = "SCAD")
  }else if(mystar=="lasso"){
    mystar.mod <-  cv.ncvreg(x, y, family = "binomial",penalty="lasso")
  }else if(mystar=="MCP"){
    mystar.mod <-  cv.ncvreg(x,y,family = "binomial", penalty = "MCP")
  }


  count =0
  count2=0
  count.local=0

  #### MCMC Iteration loop ####
  for(t in 1:(nmc+burn)){

    if (t%%1000 == 0){cat(paste(t," "))}

    #### Update Gamma and delta jointly ####
    start_time <- Sys.time()
    # Code for full enumeration when p is small
    # commented out for now; might be available in future
    # if(p<5){
    #   gam.all <- expand.grid(replicate(p, 0:1, simplify = FALSE))
    #   all.mod <- apply(gam.all[-1, ], 1, function(gam){
    #     glm(ystar~x[ ,as.logical(gam)],family = binomial(link = "logit"))})
    #   all.mod <- append(list(`1`=glm(ystar~1,family = binomial(link = "logit"))),all.mod)
    #   kap=y-n_i/2
    #   gam.logprob <- lapply(all.mod, FUN = function(mod){
    #     X.mod <- model.matrix(mod)
    #     if(exact.mixture.g==TRUE | seb.held==TRUE){
    #       mz.mod <- rep(0, length(ystar))
    #     }else{
    #       mz.mod <- mod$linear.predictors
    #     }
    #     if(seb.held==FALSE){
    #       hessian.mod <-  t(X.mod)%*%
    #         diag(mod$fitted.values*(1-mod$fitted.values))%*% X.mod
    #     }else{
    #       hessian.mod <-  0.25*(t(X.mod)%*%  X.mod)
    #     }
    #     Vz.mod <- WoodburyMatrix(A=diag(omega), B=1/delta*hessian.mod, U=X.mod,V=t(X.mod))
    #
    #
    #     z=kap/omega
    #     if(model.prior=="beta-binomial"){
    #       oj.num.log = -lchoose(p,ncol(X.mod)-1) -0.5*t(z-mz.mod)%*%
    #         solve(Vz.mod)%*% (z-mz.mod) - 0.5*
    #         determinant(Vz.mod,logarithm=TRUE)$modulus
    #     }else if (model.prior=="Uniform"){
    #       oj.num.log = -0.5*t(z-mz.mod)%*%
    #         solve(Vz.mod)%*% (z-mz.mod) - 0.5*
    #         determinant(Vz.mod,logarithm=TRUE)$modulus
    #     }
    #     return(oj.num.log)
    #
    #   })
    #   gam.logprob <- unlist(lapply(gam.logprob, as.numeric))
    #   #gam.logprob <- gam.logprob-gam.logprob[1]
    #   gam.prob <- exp(gam.logprob)/sum(exp(gam.logprob))
    #   gam <- as.matrix(gam.all[sample(1:nrow(gam.all),1,prob = gam.prob), ])
    #
    #}else{
    # Propose gamma
    gam.prop <- proposal.gamma(gam)

    # Given Gamma. propose Delta
    if(hyper==TRUE){
      # Propose delta | gam
      if(hyper.type=="robust"){
        a.prop <- (n-sum(gam.prop))/(sum(gam.prop)+1)
        a.curr <- (n-sum(gam))/(sum(gam)+1)
      }else if(hyper.type=="hyper-g"|hyper.type=="hyper-g/n"|hyper.type=="gamma"){
        a.prop <- a.curr <- 0
      }
      # Propose delta from reflective normal where a may depend on gamma
      delta.cand <- rrefnorm(1, a=a.prop, mean=delta, sd=tau)

      if(hyper.type=="gamma") {prior.d <-
        dgamma(delta.cand/n, shape = hyper.param, rate = hyper.param, log = TRUE) -
        dgamma(delta/n, shape = hyper.param, rate = hyper.param, log = TRUE)}
      if (hyper.type=="hyper-g/n") {prior.d <-
        dhypergn(delta.cand,alpha = hyper.param,n=n,logarithm = TRUE)-
        dhypergn(delta,alpha = hyper.param,n=n,logarithm = TRUE)}
      if (hyper.type=="hyper-g") {prior.d <-
        dhypergn(delta.cand,alpha = hyper.param,n=1,logarithm = TRUE)-
        dhypergn(delta,alpha = hyper.param,n=1,logarithm = TRUE)}
      if(hyper.type=="robust"){prior.d <-
        drobust(delta.cand, n=n, rho=1/(sum(gam.prop)+1),logarithm = TRUE)-
        drobust(delta, n=n, rho=1/(sum(gam)+1),logarithm = TRUE)}

    }else{
      delta.cand <- delta
    }

    x.prop = x[ ,as.logical(gam.prop)]
    x.curr = x[ ,as.logical(gam)]

    if(sum(gam.prop)==0){
      m.prop <- glm(ystar~1,family = binomial(link = "logit"))
      m.prop.data <- glm(y~1,family = binomial(link = "logit"))
    }else{
      m.prop <- glm(ystar~x.prop,family = binomial(link = "logit"))
      m.prop.data <- glm(y~x.prop,family = binomial(link = "logit"))
    }


    if(sum(gam)==0){
      m.curr <- glm(ystar~1,family = binomial(link = "logit"))
      m.curr.data <- glm(y~1,family = binomial(link = "logit"))
    }else{
      m.curr <- glm(ystar~x.curr,family = binomial(link = "logit"))
      m.curr.data <- glm(y~x.curr,family = binomial(link = "logit"))
    }

    X.prop <- model.matrix(m.prop) # with intercept
    X.curr <- model.matrix(m.curr) # with intercept
    if(exact.mixture.g==TRUE | seb.held==TRUE){
      b.prop <- rep(0, length(m.prop$coefficients))
      b.curr <- rep(0, length(m.curr$coefficients))
    }else{
      b.prop <- m.prop$coefficients
      b.curr <- m.curr$coefficients
    }
    if(seb.held==FALSE){
      hessian.prop <-  t(X.prop)%*% diag(m.prop$fitted.values*(1-m.prop$fitted.values))%*% X.prop
      hessian.curr <-  t(X.curr)%*% diag(m.curr$fitted.values*(1-m.curr$fitted.values))%*% X.curr
    }else{
      hessian.prop <-  0.25*(t(X.prop)%*%  X.prop)
      hessian.curr <-  0.25*(t(X.curr)%*% X.curr)
    }

    H.prop.data <- t(X.prop)%*% diag(m.prop.data$fitted.values*(1-m.prop.data$fitted.values))%*% X.prop
    H.curr.data <- t(X.curr)%*% diag(m.curr.data$fitted.values*(1-m.curr.data$fitted.values))%*% X.curr
    b.prop.data <- m.prop.data$coefficients
    b.curr.data <- m.curr.data$coefficients

    d.gam.prop <- H.prop.data%*% b.prop.data + 1/delta.cand * hessian.prop %*% b.prop
    d.gam.curr <- H.curr.data%*% b.curr.data + 1/delta * hessian.curr %*% b.curr

    E.gam.prop <- H.prop.data+1/delta.cand*hessian.prop
    E.gam.curr <- H.curr.data+1/delta*hessian.curr
    E.gam.prop.inv <- solve(E.gam.prop)
    E.gam.curr.inv <- solve(E.gam.curr)

    loglik.prop <- sum(y*log(m.prop.data$fitted.values)+(1-y)*log(1-m.prop.data$fitted.values))
    loglik.curr <- sum(y*log(m.curr.data$fitted.values)+(1-y)*log(1-m.curr.data$fitted.values))
    mgam.ysd.prop <-  loglik.prop +
      0.5*determinant(hessian.prop,logarithm = TRUE)$modulus-
      0.5*(sum(gam.prop)+1)*log(delta.cand)-
      0.5*determinant(E.gam.prop,logarithm = TRUE)$modulus-
      0.5*(t(b.prop.data)%*% H.prop.data%*%b.prop.data+
             1/delta.cand*t(b.prop)%*% hessian.prop%*%b.prop-
             t(d.gam.prop)%*%E.gam.prop.inv%*%d.gam.prop)
    mgam.ysd.curr <-  loglik.curr +
      0.5*determinant(hessian.curr,logarithm = TRUE)$modulus-
      0.5*(sum(gam)+1)*log(delta)-
      0.5*determinant(E.gam.curr,logarithm = TRUE)$modulus-
      0.5*(t(b.curr.data)%*% H.curr.data%*%b.curr.data+
             1/delta*t(b.curr)%*% hessian.curr%*%b.curr-
             t(d.gam.curr)%*%E.gam.curr.inv%*%d.gam.curr)

    #Vz.prop <- WoodburyMatrix(A=diag(omega), B=1/delta.cand*hessian.prop, U=X.prop,V=t(X.prop))
    #Vz.curr <- WoodburyMatrix(A=diag(omega), B=1/delta*hessian.curr, U=X.curr,V=t(X.curr))

    #kap=y-n_i/2
    #z=kap/omega

    # Calculate acceptance probability for (gam.prop,delta.cand)
    # under beta-binomial
    if(model.prior=="beta-binomial"){
      oj.num.log = -lchoose(p,sum(gam.prop))+
        mgam.ysd.prop+
        ifelse(hyper==TRUE,prior.d+
                 drefnorm(delta, a=a.curr, mean=delta.cand,sd=tau,logarithm=TRUE) ,0)

      oj.den.log = -lchoose(p,sum(gam))+
        mgam.ysd.curr+
        ifelse(hyper==TRUE,drefnorm(delta.cand, a=a.prop, mean=delta,sd=tau,logarithm=TRUE),0)
    }else if (model.prior=="Uniform"){ # Under Uniform model prior
      oj.num.log =  mgam.ysd.prop+
        ifelse(hyper==TRUE,prior.d+
                 drefnorm(delta, a=a.curr, mean=delta.cand,sd=tau,logarithm=TRUE),0)
      oj.den.log =  mgam.ysd.curr+
        ifelse(hyper==TRUE,drefnorm(delta.cand, a=a.prop, mean=delta,sd=tau,logarithm=TRUE),0)
    }
    u.gam <- runif(1)
    a.gam.prob <- min(as.numeric(exp(oj.num.log-oj.den.log)),1)
    if(u.gam<=a.gam.prob){
      gam <- gam.prop
      delta <- delta.cand
      mgam.ysd <- mgam.ysd.prop
      E.gam <- E.gam.prop
      E.gam.inv <- E.gam.prop.inv
      d.gam <- d.gam.prop
      X.inter.gam <- X.prop
    }else{
      gam <- gam
      delta <- delta
      mgam.ysd <- mgam.ysd.curr
      E.gam <- E.gam.curr
      E.gam.inv <- E.gam.curr.inv
      d.gam <- d.gam.curr
      X.inter.gam <- X.curr
    }
    #    }

    end_time <- Sys.time()
    timemat[t,1] <- end_time-start_time


    #### Update Beta #####
    if(sample.beta==TRUE){


    start_time <- Sys.time()
    # if(sum(gam)==0){
    #   mod.gam <- glm(ystar~1,family = binomial(link = "logit"))
    # }else{
    #   x.gam <- x[ , as.logical(gam)]
    #   mod.gam <- glm(ystar~x.gam,family = binomial(link = "logit"))
    # }
    #
    # X.inter.gam <- model.matrix(mod.gam) # with intercept
    # if(seb.held==TRUE){
    #   hess.mat <- 0.25*(t(X.inter.gam) %*% X.inter.gam)
    #   prior.vcov <- 4*solve(t(X.inter.gam) %*% X.inter.gam)
    # }else{
    #   hess.mat <- t(X.inter.gam) %*% diag(mod.gam$fitted.values*(1-mod.gam$fitted.values)) %*% X.inter.gam#solve(vcov(mod.gam))#
    #   prior.vcov <- vcov(mod.gam)
    # }
    #
    #
    # if(exact.mixture.g==TRUE| seb.held==TRUE){
    #   beta.hat <- rep(0, ncol(X.inter.gam))
    # }else{
    #   beta.hat <- mod.gam$coefficients
    # }
    #
    # V.omega.inv <- 1/delta*hess.mat+t(X.inter.gam)%*%diag(omega)%*% X.inter.gam
    # # Alternate version using Woodbury Matrix package
    # #V.omega.inv <- WoodburyMatrix(A=delta*prior.vcov, B=diag(omega^(-1)),U=t(X.inter.gam),V=X.inter.gam)
    # V.omega=solve(V.omega.inv)
    # V.omega <-as.matrix((V.omega+t(V.omega))/2)
    # m.omega= V.omega%*%(t(X.inter.gam) %*% kap + (1/delta)*hess.mat%*%beta.hat)
    # b = t(rmvnorm(1, mean=m.omega, sigma = V.omega))
    b=t(rmvnorm(1,mean = E.gam.inv %*% d.gam, sigma = E.gam.inv))
    end_time <- Sys.time()
    timemat[t,2] <- end_time-start_time


    }
    start_time <- Sys.time()
    # #### Update Omega ####
    # omega=unlist(mapply(rpg, n_i, X.inter.gam%*%b, num=1))
    end_time <- Sys.time()
    timemat[t,3] <- end_time-start_time

    #### Update y^* ####
    start_time <- Sys.time()
    #tic("ystar loop")
    if(exact.mixture.g==FALSE & seb.held==FALSE){

      pi0= expit(b[1])
      if(sum(gam)==0){
        pi.gam = 1/2
      }else{
        pi.gam = expit(as.matrix(X.inter.gam[ ,-1])%*% as.matrix(b[-1]))
      }

      a0 <- 1/n*log(pi0)+1/delta*log(pi.gam)
      b0 <- 1/n*log(1-pi0)+1/delta*log(1-pi.gam)
      pi.star <- expit(a0-b0)

      prop.obj <-  proposal.ystar(ystar, pi.star, n_i)
      y.star.cand <- prop.obj$y.star.cand
      local <- prop.obj$s # s=1 implies local vs s=2 implies global

      # Check if proposed ystar causes separation
      #tic("sep check")
      #sep.time[t,1] <- system.time({
      m.full <- glm(y.star.cand~x,family = binomial(link = "logit"),method = "detect_separation", solver="glpk")
      #})[3]
      #toc()

      # If yes, continue with existing one
      if(m.full$separation==TRUE){
        ystar <- ystar
      }else{
        # Else calculate acceptance probability for proposed imaginary sample
        a.prob <- ystar.acceptance.prob(y.star.cand,ystar,x,gam,b, delta, n_i, pi.star,local,mystar.mod)
        if(runif(1)<=a.prob){
          count=count+1
          if(local==1){count.local=count.local+1}
          ystar <- y.star.cand
        }else{
          ystar <- ystar
        }
      }

    }
    #toc()
    end_time <- Sys.time()
    timemat[t,4] <- end_time-start_time

    #### Update delta ####
    start_time <- Sys.time()

    if(hyper==TRUE){
      # Propose delta | current value of gamma
      # decide choice of left reflection boundary for reflective Gaussian
      # based on hyper prior type
      if(hyper.type=="robust"){
        a.curr <- (n-sum(gam))/(sum(gam)+1)
      }else if(hyper.type=="hyper-g"|hyper.type=="hyper-g/n"|hyper.type=="gamma"){
        a.prop <- a.curr <- 0
      }
      # Propose new value of delta
      delta.cand <- rrefnorm(1, a=a.prop, mean=delta, sd=tau)
      # Calculate acceptance probability for the new value of delta
      a.prob.delta <- delta.acceptance.prob(delta.cand,delta,b,ystar,x,gam,hyper.type,hyper.param,exact.mixture.g,seb.held)
      if(runif(1)<=a.prob.delta){
        count2=count2+1
        delta <- delta.cand
      }else{
        delta <- delta
      }

    }
    end_time <- Sys.time()
    timemat[t,5] <- end_time-start_time


    betastore=rep(0,p+1)

    #Save results post burn-in
    if(t > burn){
      GammaSave[t-burn, ] = gam
      count3=1
      gam.withintercept <- c(1,gam)
      for (k in 1:(p+1)){

        if(gam.withintercept[k]==1){
          betastore[k]= b[count3]
          count3=count3+1
        }
      }
      BetaSave[t-burn, ] = betastore
      #OmegaSave[t-burn, ] = omega
      ystarSave[t-burn, ] = ystar
      deltaSave[t-burn, 1] = delta
    }

  }

  # Store results as a list
  result <- list("BetaSamples"=BetaSave,
                 "GammaSamples"=GammaSave,
                 "OmegaSamples"=OmegaSave,
                 "ystarSamples"=ystarSave,
                 "deltaSamples"=deltaSave,
                 "acc.ratio.ystar"= count/(nmc+burn),
                 "acc.ratio.ystar.local"=count.local/(nmc+burn),
                 "acc.ratio.delta"=count2/(nmc+burn),
                 "timemat"=timemat)
  # Return object
  return(result)

}



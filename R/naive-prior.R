#### Load libraries #####

library(BayesLogit)
library(mvtnorm)
library(robustbase)
library(WoodburyMatrix)
library(testit)
library(detectseparation)
library(ncvreg)

#### Helper functions ####


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



naive.prior <- function(x,y,
                        prior.sdev,
                        burn=1000,
                        nmc=5000,
                        model.prior="beta-binomial",
                        exact=TRUE){

  n <- length(y)
  p <- ncol(x)


  # create matrix to save variables
  GammaSave = matrix(NA, nmc, p)
  BetaSave = matrix(NA, nmc, p+1)
  OmegaSave = matrix(NA, nmc, n)
  timemat <- matrix(NA, nmc+burn, 4)

  # Intialize parameters
  gam = rep(0,p)
  ful.mod <- glm(y~1,family = binomial(link = "logit"))
  b = c(ful.mod$coefficients,rep(0,p))
  n_i = rep(1,n)
  omega = unlist(lapply(n_i, rpg,num=1,z=0)) # prior is PG(n_i, 0) where n_i=1 for binary logit

  #### MCMC Iteration loop ####
  for(t in 1:(nmc+burn)){

    if (t%%1000 == 0){cat(paste(t," "))}

    start_time <- Sys.time()
    # Propose gamma
    gam.prop <- proposal.gamma(gam)

    x.prop = x[ ,as.logical(gam.prop)]
    x.curr = x[ ,as.logical(gam)]

    kap=y-n_i/2
    z=kap/omega

    if(sum(gam.prop)==0){
      X.prop <- model.matrix(y~1)
    }else{
      X.prop <- model.matrix(y~x.prop)
    }

    if(sum(gam)==0){
      X.curr <- model.matrix(y~1)
    }else{
      X.curr <- model.matrix(y~x.curr)
    }


    ##### Update gamma
    if(exact==TRUE){

      # Calculate acceptance probability for (gam.prop,delta.cand)
      Vz.prop <- WoodburyMatrix(A=diag(omega), B=prior.sdev^2*diag(sum(gam.prop)+1), U=X.prop,V=t(X.prop))
      Vz.curr <- WoodburyMatrix(A=diag(omega), B=prior.sdev^2*diag(sum(gam)+1), U=X.curr,V=t(X.curr))


      # under beta-binomial
      if(model.prior=="beta-binomial"){
        oj.num.log = -lchoose(p,sum(gam.prop)) -0.5*
          t(z)%*% solve(Vz.prop)%*% (z) -
          0.5*determinant(Vz.prop,logarithm=TRUE)$modulus

        oj.den.log = -lchoose(p,sum(gam)) -0.5*
          t(z)%*%solve(Vz.curr)%*% (z) -
          0.5*determinant(Vz.curr,logarithm=TRUE)$modulus
      }else if (model.prior=="Uniform"){ # Under Uniform model prior
        oj.num.log =  -0.5*t(z)%*% solve(Vz.prop)%*% (z) -
          0.5*determinant(Vz.prop,logarithm=TRUE)$modulus
        oj.den.log =  -0.5*t(z)%*%solve(Vz.curr)%*% (z) -
          0.5*determinant(Vz.curr,logarithm=TRUE)$modulus
      }



    }else{

      if(sum(gam.prop)==0){
        #m.prop <- glm(ystar~1,family = binomial(link = "logit"))
        m.prop.data <- glm(y~1,family = binomial(link = "logit"))
      }else{
        #m.prop <- glm(ystar~x.prop,family = binomial(link = "logit"))
        m.prop.data <- glm(y~x.prop,family = binomial(link = "logit"))
      }


      if(sum(gam)==0){
        #m.curr <- glm(ystar~1,family = binomial(link = "logit"))
        m.curr.data <- glm(y~1,family = binomial(link = "logit"))
      }else{
        #m.curr <- glm(ystar~x.curr,family = binomial(link = "logit"))
        m.curr.data <- glm(y~x.curr,family = binomial(link = "logit"))
      }

      loglik.prop <- sum(y*log(m.prop.data$fitted.values)+(1-y)*log(1-m.prop.data$fitted.values))
      loglik.curr <- sum(y*log(m.curr.data$fitted.values)+(1-y)*log(1-m.curr.data$fitted.values))

      H.prop.data <- t(X.prop)%*% diag(m.prop.data$fitted.values*(1-m.prop.data$fitted.values))%*% X.prop
      H.curr.data <- t(X.curr)%*% diag(m.curr.data$fitted.values*(1-m.curr.data$fitted.values))%*% X.curr
      b.prop.data <- m.prop.data$coefficients
      b.curr.data <- m.curr.data$coefficients

      d.gam.prop <- H.prop.data%*% b.prop.data
      d.gam.curr <- H.curr.data%*% b.curr.data

      E.gam.prop <- H.prop.data+ 1/prior.sdev^2*diag(sum(gam.prop)+1)
      E.gam.curr <- H.curr.data+1/prior.sdev^2*diag(sum(gam)+1)
      E.gam.prop.inv <- solve(E.gam.prop)
      E.gam.curr.inv <- solve(E.gam.curr)

      mgam.ysd.prop <-  loglik.prop -
        (sum(gam.prop)+1)*log(prior.sdev)-
        0.5*determinant(E.gam.prop,logarithm = TRUE)$modulus-
        0.5*(t(b.prop.data)%*% H.prop.data%*%b.prop.data-
               t(d.gam.prop)%*%E.gam.prop.inv%*%d.gam.prop)
      mgam.ysd.curr <-  loglik.curr -
        0.5*(sum(gam)+1)*log(prior.sdev)-
        0.5*determinant(E.gam.curr,logarithm = TRUE)$modulus-
        0.5*(t(b.curr.data)%*% H.curr.data%*%b.curr.data-
               t(d.gam.curr)%*%E.gam.curr.inv%*%d.gam.curr)



      # Calculate acceptance probability for (gam.prop,delta.cand)
      # under beta-binomial
      if(model.prior=="beta-binomial"){
        oj.num.log = -lchoose(p,sum(gam.prop))+
          mgam.ysd.prop

        oj.den.log = -lchoose(p,sum(gam))+
          mgam.ysd.curr

      }else if (model.prior=="Uniform"){ # Under Uniform model prior
        oj.num.log =  mgam.ysd.prop
        oj.den.log =  mgam.ysd.curr
      }


    }

    u.gam <- runif(1)
    a.gam.prob <- min(as.numeric(exp(oj.num.log-oj.den.log)),1)
    if(u.gam<=a.gam.prob){
      gam <- gam.prop
      X.inter.gam <- X.prop
      if(exact==FALSE){
        E.gam <- E.gam.prop
        E.gam.inv <- E.gam.prop.inv
        d.gam <- d.gam.prop
      }


      #delta <- delta.cand
    }else{
      gam <- gam
      X.inter.gam <- X.curr
      if(exact==FALSE){
        E.gam <- E.gam.curr
        E.gam.inv <- E.gam.curr.inv
        d.gam <- d.gam.curr

      }

      #delta <- delta
    }
    #    }
    end_time <- Sys.time()
    timemat[t,3] <- end_time-start_time

    #### Update Beta #####
    start_time <- Sys.time()

    #### Update Beta and omega (if exact=TRUE)
    if(exact==TRUE){
      V.omega.inv <- 1/prior.sdev^2*diag(sum(gam)+1)+t(X.inter.gam)%*%diag(omega)%*% X.inter.gam
      # Alternate version using Woodbury Matrix package
      #V.omega.inv <- WoodburyMatrix(A=delta*prior.vcov, B=diag(omega^(-1)),U=t(X.inter.gam),V=X.inter.gam)
      V.omega=solve(V.omega.inv)
      V.omega <-as.matrix((V.omega+t(V.omega))/2)
      m.omega= V.omega%*%(t(X.inter.gam) %*% kap)

      b = t(rmvnorm(1, mean=m.omega, sigma = V.omega))



      #### Update Omega ####
      start_time <- Sys.time()
      omega=unlist(mapply(rpg, n_i, X.inter.gam%*%b, num=1))
      end_time <- Sys.time()
      timemat[t,3] <- end_time-start_time


    }else{


      b=t(rmvnorm(1,mean = E.gam.inv %*% d.gam, sigma = E.gam.inv))


    }
    end_time <- Sys.time()
    timemat[t,2] <- end_time-start_time


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
      if(exact==TRUE){
        OmegaSave[t-burn, ] = omega
      }

      #ystarSave[t-burn, ] = ystar
      #deltaSave[t-burn, 1] = delta
    }


  }

  # Store results as a list
  result <- list("BetaSamples"=BetaSave,
                 "GammaSamples"=GammaSave,
                 "OmegaSamples"=OmegaSave,
                 #"ystarSamples"=ystarSave,
                 #"deltaSamples"=deltaSave,
                 #"acc.ratio.ystar"= count/(nmc+burn),
                 #"acc.ratio.ystar.local"=count.local/(nmc+burn),
                 #"acc.ratio.delta"=count2/(nmc+burn),
                 "timemat"=timemat)
  # Return object
  return(result)



}

library(BayesLogit)
library(mvtnorm)
library(robustbase)
library(WoodburyMatrix)
library(testit)
library(detectseparation)



rrefnorm <- function(n,a=0, mean=0, sd=1){
  eps <- rnorm(n)
  g=a+abs(mean+sd*eps-a)
  return(g)
}

drefnorm <- function(x,a=0,mean=0,sd=1,logarithm=FALSE){
  if(x>=a){
    val <- dnorm(x,mean=mean,sd=sd)+dnorm(x,mean=2*a-mean, sd=sd)
  }else{
    val <- 0
  }
  if(logarithm==TRUE){return(log(val))}
  return(val)
}


dhypergn <- function(g, alpha=4,n,logarithm=FALSE){
  val <- (alpha-2)/(2*n)*(1/(1+g/n))^(alpha/2)
  if(logarithm==TRUE){
    val <- log(alpha-2)-log(2*n)-alpha/2*log(1+g/n)
  }
  return (val)
}


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



ystar.acceptance.prob <- function(y.new,y.old,xmat,gam, bmat,delta,n_i,pi.star,local){
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


  new.log.bgam.dens <- dmvnorm(t(bmat), mean = mod.new$coefficients, sigma = delta*vcov(mod.new),log=TRUE)
  old.log.bgam.dens <- dmvnorm(t(bmat), mean = mod.old$coefficients, sigma = delta*vcov(mod.old),log=TRUE)

  log.num <- lgamma(new.sum+0.5)+lgamma(n-new.sum+0.5) + new.log.bgam.dens
  +ifelse(local==2,dbinom(y.old,n_i,prob=pi.star,log=TRUE),0)#+log(dbinom(y.old,1,prob=0.5))
  log.den <- lgamma(old.sum+0.5)+lgamma(n-old.sum+0.5) + old.log.bgam.dens
  +ifelse(local==2,dbinom(y.new,n_i,prob=pi.star,log=TRUE),0)#+log(dbinom(y.new,1,prob=0.5))

  a.p <- min(exp(log.num-log.den),1)

  return (a.p)

}

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
    #+ dlnorm(delta.old, meanlog = log(delta.new),sdlog = 0.25, log = TRUE)#log(dgamma(delta.old,delta.new,1))
  log.den <- dmvnorm(t(bmat), mean = beta.star, sigma = delta.old*inf.mat,log=TRUE)
    #+ dlnorm(delta.new, meanlog = log(delta.old),sdlog = 0.25, log = TRUE)#log(dgamma(delta.new,delta.old,1))

  a.p.delta <- min(exp(log.num-log.den),1)

  return (a.p.delta)
}

resample <- function(x, ...) x[sample.int(length(x), ...)]

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



  }else if(s==2){ # swap one entry from active set with one entry from non-active set randomly
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


proposal.ystar <- function(y.star, pi.star, n_i, d=5){
  n <- length(y.star)
  prob1 <- c(0.5,0.2,0.15,0.10,0.05)
  prob2 <- c(0.7,0.3)

  s <- resample(1:2,1, prob = prob2)
  # Use technique of George mchlloch stat sinaca paper  eq 46 if s==1
  if(s==1){
    d0 <- resample(1:d,1, prob=prob1)
    ind <- resample(1:n,d0)
    y.star.cand <- y.star
    y.star.cand[ind]=as.numeric(!y.star.cand[ind])

  }else if(s==2){
    
    if(any(is.nan(pi.star))){
      s=1
      d0 <- resample(1:d,1, prob=prob1)
      ind <- resample(1:n,d0)
      y.star.cand <- y.star
      y.star.cand[ind]=as.numeric(!y.star.cand[ind])
    }else{
      y.star.cand <- rbinom(n, n_i,  pi.star)  
    }
    
  }
  return(list(y.star.cand=y.star.cand,s=s))
}



expit <- function(x){
  return(1/(1+exp(-x)))
}

rhypergn <- function(nsamp, alpha=4,n){
  u <- rbeta(nsamp, alpha/2-1,1)
  return (n*u/(1-u))
}

Laplace.pep <- function(x,y,burn=1000,nmc=5000, model.prior="beta-binomial",hyper="TRUE", hyper.type="gamma",hyper.param=NULL,exact.mixture.g=FALSE,seb.held=FALSE){

  n <- length(y)
  p <- ncol(x)
  tau=n/2
  # create matrix to save variables
  GammaSave = matrix(NA, nmc, p)
  BetaSave = matrix(NA, nmc, p+1)
  OmegaSave = matrix(NA, nmc, n)
  ystarSave = matrix(NA, nmc, n)
  deltaSave = matrix(NA,nmc, 1)
  timemat <- matrix(NA, nmc+burn, 5)

  # Intialize
  gam = rep(0,p)
  ful.mod <- glm(y~1,family = binomial(link = "logit"))
  b = c(ful.mod$coefficients,rep(0,p))
  n_i = rep(1,n)
  omega = unlist(lapply(n_i, rpg,num=1,z=0)) # prior is PG(n_i, 0) where n_i=1 for binary logit

  if(exact.mixture.g==TRUE){
    ystar=y
  }else{
    ystar = rbinom(n,1,0.5)
  }
  delta = n#rhypergn(1,alpha = 4,n=n)


  count =0
  count2=0
  count.local=0

  for(t in 1:(nmc+burn)){

    if (t%%1000 == 0){cat(paste(t," "))}

    #### Update Gamma and delta jointly ####
    start_time <- Sys.time()
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

      if(hyper==TRUE){
        # Propose delta | gam
        if(hyper.type=="robust"){
          a.prop <- (n-sum(gam.prop))/(sum(gam.prop)+1)
          a.curr <- (n-sum(gam))/(sum(gam)+1)
        }else if(hyper.type=="hyper-g"|hyper.type=="hyper-g/n"|hyper.type=="gamma"){
          a.prop <- a.curr <- 0
        }

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
      }else{
        m.prop <- glm(ystar~x.prop,family = binomial(link = "logit"))
      }


      if(sum(gam)==0){
        m.curr <- glm(ystar~1,family = binomial(link = "logit"))
      }else{
        m.curr <- glm(ystar~x.curr,family = binomial(link = "logit"))
      }

      X.prop <- model.matrix(m.prop) # with intercept
      X.curr <- model.matrix(m.curr) # with intercept
      if(exact.mixture.g==TRUE | seb.held==TRUE){
        mz.prop <- rep(0, length(y))
        mz.curr <- rep(0, length(y))
      }else{
        mz.prop <- m.prop$linear.predictors
        mz.curr <- m.curr$linear.predictors
      }
      if(seb.held==FALSE){
        hessian.prop <-  t(X.prop)%*% diag(m.prop$fitted.values*(1-m.prop$fitted.values))%*% X.prop
        hessian.curr <-  t(X.curr)%*% diag(m.curr$fitted.values*(1-m.curr$fitted.values))%*% X.curr
      }else{
        hessian.prop <-  0.25*(t(X.prop)%*%  X.prop)
        hessian.curr <-  0.25*(t(X.curr)%*% X.curr)
      }

      Vz.prop <- WoodburyMatrix(A=diag(omega), B=1/delta.cand*hessian.prop, U=X.prop,V=t(X.prop))
      Vz.curr <- WoodburyMatrix(A=diag(omega), B=1/delta*hessian.curr, U=X.curr,V=t(X.curr))
      #Vz.prop <- diag(omega^{-1})+delta.cand*(X.prop %*% vcov(m.prop) %*% t(X.prop))
      #Vz.curr <- diag(omega^{-1})+delta*(X.curr %*% vcov(m.curr) %*% t(X.curr))

      kap=y-n_i/2
      z=kap/omega


      if(model.prior=="beta-binomial"){
        oj.num.log = -lchoose(p,sum(gam.prop)) -0.5*
          t(z-mz.prop)%*% solve(Vz.prop)%*% (z-mz.prop) -
          0.5*determinant(Vz.prop,logarithm=TRUE)$modulus+
          ifelse(hyper==TRUE,prior.d+
                   drefnorm(delta, a=a.curr, mean=delta.cand,sd=tau,logarithm=TRUE) ,0)

        oj.den.log = -lchoose(p,sum(gam)) -0.5*
          t(z-mz.curr)%*%solve(Vz.curr)%*% (z-mz.curr) -
          0.5*determinant(Vz.curr,logarithm=TRUE)$modulus+
          ifelse(hyper==TRUE,drefnorm(delta.cand, a=a.prop, mean=delta,sd=tau,logarithm=TRUE),0)
      }else if (model.prior=="Uniform"){
        oj.num.log =  -0.5*t(z-mz.prop)%*% solve(Vz.prop)%*% (z-mz.prop) -
          0.5*determinant(Vz.prop,logarithm=TRUE)$modulus+
          ifelse(hyper==TRUE,prior.d+
                   drefnorm(delta, a=a.curr, mean=delta.cand,sd=tau,logarithm=TRUE),0)
        oj.den.log =  -0.5*t(z-mz.curr)%*%solve(Vz.curr)%*% (z-mz.curr) -
          0.5*determinant(Vz.curr,logarithm=TRUE)$modulus+
          ifelse(hyper==TRUE,drefnorm(delta.cand, a=a.prop, mean=delta,sd=tau,logarithm=TRUE),0)
      }
      #oj=as.numeric(exp(oj.num.log-oj.den.log))
      #print(oj)
      u.gam <- runif(1)
      a.gam.prob <- min(as.numeric(exp(oj.num.log-oj.den.log)),1)
      if(u.gam<=a.gam.prob){
        gam <- gam.prop
        delta <- delta.cand
      }else{
        gam <- gam
        delta <- delta
      }
#    }

    end_time <- Sys.time()
    timemat[t,1] <- end_time-start_time


    #### Update Beta #####
    start_time <- Sys.time()
    if(sum(gam)==0){
      mod.gam <- glm(ystar~1,family = binomial(link = "logit"))
    }else{
      x.gam <- x[ , as.logical(gam)]
      mod.gam <- glm(ystar~x.gam,family = binomial(link = "logit"))
    }

    X.inter.gam <- model.matrix(mod.gam) # with intercept
    if(seb.held==TRUE){
      hess.mat <- 0.25*(t(X.inter.gam) %*% X.inter.gam)
      prior.vcov <- 4*solve(t(X.inter.gam) %*% X.inter.gam)
    }else{
      hess.mat <- t(X.inter.gam) %*% diag(mod.gam$fitted.values*(1-mod.gam$fitted.values)) %*% X.inter.gam#solve(vcov(mod.gam))#
      prior.vcov <- vcov(mod.gam)
    }


    if(exact.mixture.g==TRUE| seb.held==TRUE){
      beta.hat <- rep(0, ncol(X.inter.gam))
    }else{
      beta.hat <- mod.gam$coefficients
    }

    V.omega.inv <- 1/delta*hess.mat+t(X.inter.gam)%*%diag(omega)%*% X.inter.gam
    #V.omega.inv <- WoodburyMatrix(A=delta*prior.vcov, B=diag(omega^(-1)),U=t(X.inter.gam),V=X.inter.gam)
    V.omega=solve(V.omega.inv)
    V.omega <-as.matrix((V.omega+t(V.omega))/2)
    m.omega= V.omega%*%(t(X.inter.gam) %*% kap + (1/delta)*hess.mat%*%beta.hat)
    b = t(rmvnorm(1, mean=m.omega, sigma = V.omega))

    end_time <- Sys.time()
    timemat[t,2] <- end_time-start_time

    start_time <- Sys.time()

    #### Update Omega ####
    omega=unlist(mapply(rpg, n_i, X.inter.gam%*%b, num=1))
    end_time <- Sys.time()
    timemat[t,3] <- end_time-start_time

    #### Update y^* ####
    start_time <- Sys.time()
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

      m.full <- glm(y.star.cand~x,family = binomial(link = "logit"),method = "detect_separation", solver="glpk")

      if(m.full$separation==TRUE){
        #print("Hello")
        ystar <- ystar
        #warncount=FALSE
      }else{
        #ystar.acceptance.prob <- function(y.new,y.old,xmat, bmat,delta,n_i,pi.star)
        a.prob <- ystar.acceptance.prob(y.star.cand,ystar,x,gam,b, delta, n_i, pi.star,local)
        if(runif(1)<=a.prob){
          count=count+1
          if(local==1){count.local=count.local+1}
          ystar <- y.star.cand
        }else{
          ystar <- ystar
        }
      }
    }

    end_time <- Sys.time()
    timemat[t,4] <- end_time-start_time

    #### Update delta ####
    start_time <- Sys.time()

    if(hyper==TRUE){
      # Propose delta | gam
      if(hyper.type=="robust"){
        a.curr <- (n-sum(gam))/(sum(gam)+1)
      }else if(hyper.type=="hyper-g"|hyper.type=="hyper-g/n"|hyper.type=="gamma"){
        a.prop <- a.curr <- 0
      }

      delta.cand <- rrefnorm(1, a=a.prop, mean=delta, sd=tau)

      #delta.acceptance.prob <- function(delta.new,delta.old,bmat, y.star,xmat)
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

    #Save results
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
      OmegaSave[t-burn, ] = omega
      ystarSave[t-burn, ] = ystar
      deltaSave[t-burn, 1] = delta
    }

  }

  result <- list("BetaSamples"=BetaSave,
                 "GammaSamples"=GammaSave,
                 "OmegaSamples"=OmegaSave,
                 "ystarSamples"=ystarSave,
                 "deltaSamples"=deltaSave,
                 "acc.ratio.ystar"= count/(nmc+burn),
                 "acc.ratio.ystar.local"=count.local/(nmc+burn),
                 "acc.ratio.delta"=count2/(nmc+burn),
                 "timemat"=timemat)

  return(result)

}



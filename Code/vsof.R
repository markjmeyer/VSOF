#### Dependent libraries ####
library(splines)
library(Matrix)
library(numbers)

#### KL lower bound for Variational Additive Models ####
eval_logp_gauss  <- function(Sigmaq, sigB, Muq, p, Sigqb, Ae, Be, n, Bqe, Af, Bf, Kfi, Bqf){
  if(n/2 < 150){
    out   <- (1/2)*(p + sum(Kfi[-1])) - (n/2)*log(2*pi) - (p/2)*log(sigB) +
      (1/2)*log(det(Sigmaq)) - (1/(2*sigB))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
      Ae*log(Be) - (Ae + n/2)*log(Bqe) + log(gamma(Ae + n/2)) - log(gamma(Ae)) +
      sum(Af*log(Bf) - (Af + Kfi[-1]/2)*log(Bqf) + log(gamma(Af + Kfi[-1]/2)) - log(gamma(Af)))
  } else{
    out   <- (1/2)*(p + sum(Kfi[-1])) - (n/2)*log(2*pi) - (p/2)*log(sigB) +
      (1/2)*log(det(Sigmaq)) - (1/(2*sigB))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
      Ae*log(Be) - (Ae + n/2)*log(Bqe) - log(gamma(Ae)) +
      sum(Af*log(Bf) - (Af + Kfi[-1]/2)*log(Bqf) + log(gamma(Af + Kfi[-1]/2)) - log(gamma(Af)))
  }
  
  return(as.numeric(out))
}

eval_logp_binpr <- function(n, Y, XZ, Muq, p, Sigqb, Sigmaq, Af, Bf, Kfi, Bqf){
  ones  <- rep(1, n)
  out   <- t(Y)%*%log(pnorm(XZ%*%Muq) + 1) +
    t(ones - Y)%*%log(ones - pnorm(XZ%*%Muq) + 1) -
    (1/(2))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
    (1/2)*log(det(Sigmaq)) +
    sum(Af*log(Bf) - (Af + Kfi[-1]/2)*log(Bqf) + log(gamma(Af + Kfi[-1]/2)) - log(gamma(Af)))
}

#### construct. functions ####
construct_P  <- function(D, xi, Theta){
  diff0   <- diag(1, D, D)
  diff2   <- matrix(rep(c(1, -2, 1, rep(0, D - 2)), D - 2)[1:((D - 2) * D)], D - 2, D, 
                    byrow = TRUE)
  P0      <- t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2      <- t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat   <- xi * P0 + (1 - xi) * P2
  
  return(P.mat)
}

construct_fe  <- function(fe, ps.spline, data, Ts){
  mfe     <- as.formula(paste(deparse(fe),"- 1"))
  Zl      <- model.matrix(mfe, data = data)
  dfs     <- ps.spline$df
  degs    <- ps.spline$degree
  xi      <- ps.spline$xi
  Ks      <- vector(length = length(Ts))
  Kf      <- length(Ts)*dfs
  Th      <- O <-   vector('list', length = length(Ts))
  Zf      <- matrix(NA, nrow(Zl), Kf)
  Tsi     <- c(0, unlist(Ts))
  Kfi     <- c(0, rep(dfs, length(Ts)))
  Kcs     <- cumsum(Kfi)
  for(i in 1:length(Ts)){
    Ks[i]     <- Ts[[i]]
    x         <- 1:Ks[i]
    Th[[i]]   <- bs(x, df = dfs, intercept = FALSE, degree = degs)
    O[[i]]    <- construct_P(Ks[i], xi, Th[[i]])
    
    Zf[,(Kcs[i]+1):Kcs[i+1]]      <- Zl[,(Tsi[i] + 1):(sum(Tsi[i:(i+1)]))]%*%Th[[i]]
  }
  
  out     <- list(Zf = Zf, Th = Th, O = O)
  
  return(out)
}

construct_se  <- function(se, se.spline, data){
  mse     <- as.formula(paste(deparse(se),"- 1"))
  Zl      <- model.matrix(mse, data = data)
  dfs     <- se.spline$df
  Kf      <- ncol(Zl)*dfs
  degs    <- se.spline$degree
  slength  <- se.spline$slength
  Zf      <- matrix(NA, nrow(Zl), Kf)
  Kfi     <- c(0, rep(dfs, ncol(Zl)))
  Kcs     <- cumsum(Kfi)
  Th      <- O <- sRange  <-  vector('list', length = ncol(Zl))
  for(i in 1:ncol(Zl)){
    Zf[,(Kcs[i]+1):Kcs[i+1]] <- Zfi <- bs(as.numeric(Zl[,i]), df = dfs, intercept = FALSE, degree = degs)
    Ks      <- ncol(Zfi)
    Db      <- diff(diag(Ks))
    
    O[[i]]        <- t(Db)%*%Db
    sRange[[i]]   <- seq(range(Zl[,i])[1], range(Zl[,i])[2], length = slength)
    Th[[i]]       <- bs(sRange[[i]], df = dfs, intercept = FALSE, degree = degs)
  }
  
  out     <- list(Zf = Zf, Th = Th, O = O, sRange = sRange)
  
  return(out)
}

#### Variational Additive Models ####
vsof   <- function(fixed, se = NULL, fe = NULL, Ts = NULL, data, 
                  se.spline = NULL, fe.spline = NULL, family = gaussian(link = "identity"),
                  priors = list(Ae = 0.01, Be = 0.01, Af = c(0.01, 0.01), # edit defaults
                                Bf = c(0.01, 0.01), sigB = 1000),
                  tol = 1e-5, maxIter = 500, up = 100, dots = up/10,
                  verbose = TRUE){
  
  ## extract data, build matrices and vectors
  mf      <- model.frame(formula = fixed, data = data)
  Y       <- as.numeric(model.response(mf))
  Xc      <- model.matrix(attr(mf,"terms"), data = mf)
  if(is.null(fe.spline)){
    fe.spline   <- list(df = 8, degree = 3, xi = 0.01)
  }
  if(is.null(se.spline)){
    se.spline   <- list(df = 4, degree = 3, slength = 100)
  }
  if(!is.null(fe)){
    cfe     <- construct_fe(fe, fe.spline, data, Ts)
    
    O           <- fe.O   <- cfe$O
    Th          <- fe.Th  <- cfe$Th
    Zf          <- fe.Zf  <- cfe$Zf
    data.names  <- fe.dn  <- names(model.frame(formula = fe, data = data))
    Kfi         <- c(0, rep(fe.spline$df, length(O)))
    Kcs         <- cumsum(Kfi)
    test.flag   <- rep('fe', length(fe.O))
  }
  if(!is.null(se)){
    cse     <- construct_se(se, se.spline, data)
    
    O                 <- se.O   <- cse$O
    Th                <- se.Th  <- cse$Th
    Zf                <- se.Zf  <- cse$Zf
    se.spline$sRange  <- cse$sRange
    data.names        <- se.dn  <- names(model.frame(formula = se, data = data))
    Kfi               <- c(0, rep(se.spline$df, length(O)))
    Kcs               <- cumsum(Kfi)
    test.flag         <- rep('se', length(se.O))
  }
  if(!is.null(fe) & !is.null(se)){
    O           <- c(se.O, fe.O)
    Th          <- c(se.Th, fe.Th)
    Zf          <- cbind(se.Zf, fe.Zf)
    data.names  <- c(se.dn, fe.dn)
    Kfi         <- c(0, rep(se.spline$df, length(se.O)), rep(fe.spline$df, length(fe.O)))
    Kcs         <- cumsum(Kfi)
    test.flag   <- c(rep('se', length(se.O)), rep('fe', length(fe.O)))
  }
  n       <- length(Y) # number of subjects
  
  ### build design matrix ###
  XZ      <- cbind(Xc, Zf)
  
  ## prior values ##
  Ae    <- priors$Ae
  Be    <- priors$Be
  if(length(priors$Af) == length(O)){
    Af    <- priors$Af # edit defaults
  } else{
    Af    <- rep(0.01, length(O))
  }
  if(length(priors$Bf) == length(O)){
    Bf    <- priors$Bf # edit defaults
  } else{
    Bf    <- rep(0.01, length(O))
  }
  sigB  <- priors$sigB
  
  ## initialize algorithm ##
  Bqe         <- 1
  Bqf         <- rep(1, length = length(O))
  logp_delta	<- 0.5
  iter		    <- 0
  logp_prev	  <- 0
  
  ## set matrices ##
  XtX         <- t(Xc)%*%Xc
  CtC         <- t(XZ)%*%XZ
  p           <- ncol(Xc)
  Kf          <- ncol(Zf)
  MuMat       <- matrix(0, nrow = maxIter, ncol = p + Kf)
  SigmaMat    <- array(0, dim = c(p + Kf, p + Kf, maxIter))
  BqfMat      <- matrix(0, nrow = maxIter, ncol = length(O))
  lDelta      <- rep(logp_delta, maxIter)
  logpi       <- rep(0, maxIter)
  logpp       <- rep(0, maxIter)
  Sigq        <- vector('list', length = length(O) + 1)
  
  if(family$family == "gaussian"){
    CtY         <- t(XZ)%*%Y
    BqeVec      <- rep(0, maxIter)
    
    ## algorithm ##
    while(logp_delta > tol){
      iter	<- iter + 1
      
      ## update \Sigma_q ##
      Sigqb   <- Sigq[[1]]   <- (1/sigB)*diag(p)
      for(f in 1:length(O)){
        Sigq[[f+1]] <- ((Af[f] + (1/2)*(Kfi[f + 1]))/Bqf[f])*O[[f]]
      }
      G       <- bdiag(Sigq)
      Precq   <- as.matrix(((Ae + n/2)/Bqe)*CtC + G)
      Sigmaq  <- solve(Precq)
      
      ## update \mu_q ##
      Muq     <- (Ae + n/2)/(Bqe)*Sigmaq%*%CtY
      Ctf     <- (Ae + n/2)/(Bqe)*Sigmaq%*%t(XZ)
      
      ## update B_qe ##
      SSY     <- as.numeric(t(Y - XZ%*%Muq)%*%(Y - XZ%*%Muq))
      Bqe     <- Be + (1/2)*(SSY + sum(diag(CtC%*%Sigmaq)))
      
      ## update B_qf ##
      for(f in 1:length(O)){
        inds    <- (Kcs[f]+1):Kcs[f+1] + 1
        MltMf   <- as.numeric(t(Muq[inds])%*%Muq[inds])
        Bqf[f]  <- Bf[1] + (1/2)*(MltMf + sum(diag(Sigq[[f+1]])))
      }
      
      ## check criteria ##
      logp_i      <- eval_logp_gauss(Sigmaq = Sigmaq, sigB = sigB, Muq = Muq, p = p, 
                                     Sigqb = Sigqb, Ae = Ae, Be = Be, n = n, Bqe = Bqe, 
                                     Af = Af, Bf = Bf, Kfi = Kfi, Bqf = Bqf)
      
      logp_delta	<- abs(logp_i - logp_prev)
      
      ## track convergence ##  
      lDelta[iter]  <- logp_delta
      logpi[iter]   <- logp_i
      logpp[iter]   <- logp_prev
      logp_prev	    <- logp_i
      
      ## store estimates ##
      MuMat[iter,]      <- as.numeric(Muq)
      SigmaMat[,,iter]  <- as.matrix(Sigmaq)
      BqeVec[iter]      <- Bqe
      BqfMat[iter,]     <- Bqf
      
      ## console update ##
      if(verbose){
        if(mod(iter, dots) == 0){
          cat('.')
        }
        
        if(mod(iter, up) == 0){
          cat(paste("\n",iter,"samples completed\n"))
        }
      }
      
      if(iter == maxIter){
        warning('Max number of iterations reached; increase max'); break
      }
      
    }
  } else{
    Muqa      <- rep(0, n)
    MuqaMat   <- matrix(0, nrow = maxIter, ncol = n)
    
    ## algorithm ##
    while(logp_delta > tol){
      iter	<- iter + 1
      
      ## update CtY
      CtY <- t(XZ)%*%Muqa
      
      ## update \Sigma_q ##
      Sigqb   <- Sigq[[1]]   <- (1/sigB)*diag(p)
      for(f in 1:length(O)){
        Sigq[[f+1]] <- ((Af[f] + (1/2)*(Kfi[f + 1]))/Bqf[f])*O[[f]]
      }
      G       <- bdiag(Sigq)
      Precq   <- as.matrix(CtC + G)
      Sigmaq  <- solve(Precq)
      
      ## update \mu_q ##
      Muq     <- Sigmaq%*%CtY
      Ctf     <- Sigmaq%*%t(XZ)
      
      ## update B_qf ##
      for(f in 1:length(O)){
        inds    <- (Kcs[f]+1):Kcs[f+1] + 1
        MltMf   <- as.numeric(t(Muq[inds])%*%Muq[inds])
        Bqf[f]  <- Bf[1] + (1/2)*(MltMf + sum(diag(Sigq[[f+1]])))
      }
      
      ## update Muqa ##
      XZMuqa  <- XZ%*%Muq
      Muqa    <- XZMuqa + dnorm(XZMuqa)/((pnorm(XZMuqa)^Y)*((pnorm(XZMuqa) - 1)^(1-Y)))
      
      ## check criteria ##
      logp_i      <- eval_logp_binpr(n = n, Y = Y, XZ = XZ, Muq = Muq, p = p, 
                                     Sigqb = Sigmaq[1:p,1:p], Sigmaq = Sigmaq, 
                                     Af = Af, Bf = Bf, Kfi = Kfi, Bqf = Bqf)
      
      logp_delta	<- abs(logp_i - logp_prev)
      
      ## track convergence ##  
      lDelta[iter]  <- logp_delta
      logpi[iter]   <- logp_i
      logpp[iter]   <- logp_prev
      logp_prev	    <- logp_i
      
      ## store estimates ##
      MuMat[iter,]      <- as.numeric(Muq)
      MuqaMat[iter,]    <- as.numeric(Muqa)
      SigmaMat[,,iter]  <- as.matrix(Sigmaq)
      BqfMat[iter,]     <- Bqf
      
      ## console update ##
      if(verbose){
        if(mod(iter, dots) == 0){
          cat('.')
        }
        
        if(mod(iter, up) == 0){
          cat(paste("\n",iter,"samples completed\n"))
        }
      }
      
      if(iter == maxIter){
        warning('Max number of iterations reached; increase max'); break
      }
      
    }   
  }
  
  ## algo outs ##
  algOut  <- list(maxIter = maxIter, tol = tol, totalIter = iter, 
                  logpd = list(lDelta = lDelta, logp = logpi, prev_logp = logpp))
  
  ## outputs ##
  xcoef     <- matrix(MuMat[iter,1:p], p, 1)
  colnames(xcoef) <- "Estimate"
  rownames(xcoef) <- colnames(Xc)
  
  dmat      <- list(XZ = XZ, Xc = Xc, Zf = Zf, Ts = Ts)
  dn        <- vector('list', length = length(Th))
  for(j in 1:length(Th)){
    dn[[j]]   <- paste(c(names(mf)[1], data.names[j]), collapse = ' by ')
  }
  if(!is.null(fe)){
    fee       <- vector('list', length = length(fe.Th))
    for(j in 1:length(fe.Th)){
      inds      <- (Kcs[j]+1):Kcs[j+1] + 1
      fee[[j]]  <- as.numeric(Th[[j]]%*%MuMat[iter,inds])
    }
  } else{
    fee   <- NULL
  }
  if(!is.null(se)){
    see       <- vector('list', length = length(se.Th))
    for(j in 1:length(se.Th)){
      inds      <- (Kcs[j]+1):Kcs[j+1] + 1
      see[[j]]  <- as.numeric(Th[[j]]%*%MuMat[iter,inds])
    }
  } else{
    see   <- NULL
  }
  if(family$family == "gaussian"){
    BqeOut    <- BqeVec[iter]
    MuqaOut   <- NULL
  } else{
    BqeOut    <- NULL
    MuqaOut   <- MuqaMat[iter,]
  }
  out       <- list(VBEstimates = MuMat[iter,],
                    Muqa = MuqaOut,
                    fixedEffects = xcoef,
                    smoothEffect = see,
                    funcEffect = fee,
                    Ctf = Ctf, data.name = dn,
                    varComps = list(Sigmaq = SigmaMat[,,iter], Bqe = BqeOut, 
                                    Bqf = BqfMat[iter,]),
                    model = list(fixed = fixed, fe = fe, se = se, 
                                 data = data, p = p, Kf = Kf, Kfi = Kfi, 
                                 family = family, test.flag = test.flag,
                                 dmat = dmat, Th = Th, O = O, Y = Y),
                    priors = priors, nfun = length(test.flag),
                    algOut = algOut,
                    se.spline = se.spline,
                    fe.spline = fe.spline)
  
  
  class(out) <- 'vsof'
  return(out)
  
}


#### ZLS Test ####
zls_test <- function(object, ...){
  UseMethod("zls_test")
}

zls_test.vsof <- function(model, test = NULL){
  
  family      <- model$model$family$family
  Th          <- model$model$Th
  SigmaMat    <- model$varComps$Sigmaq
  p           <- model$model$p
  Kf          <- model$model$Kf
  Kfi         <- model$model$Kfi
  Kcs         <- cumsum(Kfi)
  Ctf         <- model$Ctf
  BqeVec      <- model$varComps$Bqe
  Muqb        <- model$VBEstimates[c(1:model$model$p)]
  Muqa        <- model$Muqa
  Ae          <- model$priors$Ae
  Y           <- model$model$Y
  test.flag   <- model$model$test.flag
  n           <- length(Y)
  BqfVec      <- model$varComps$Bqf
  XZ          <- model$model$dmat$Xc
  out         <- vector('list', length = length(Th))
  for(f in 1:length(Th)){
    
    ## ZLS test ##
    Sigmaqf   <- SigmaMat[(Kcs[f]+1):Kcs[f+1],(Kcs[f]+1):Kcs[f+1]]
    ctf       <- Th[[f]]%*%Ctf[(Kcs[f]+1):Kcs[f+1],]
    
    if(family == "gaussian"){
      V0        <- ((BqeVec[1])/(Ae + n/2 - 1))*diag(n)
    } else{
      XZMuqa    <- XZ%*%Muqb
      sig0      <- sqrt(1 + sum(BqfVec))
      Vau0      <- 1 + ((XZMuqa/sig0)*dnorm(-(XZMuqa/sig0))/pnorm(-(XZMuqa/sig0))) - (dnorm(-(XZMuqa/sig0))/pnorm(-(XZMuqa/sig0)))^2
      Vau1      <- 1 - ((XZMuqa/sig0)*dnorm(-(XZMuqa/sig0))/(1 - pnorm(-(XZMuqa/sig0)))) - (dnorm(-(XZMuqa/sig0))/(1 - pnorm(-(XZMuqa/sig0))))^2
      Vars      <- ((1*(sign(Muqa) < 0))*Vau0 + (1*(sign(Muqa) >= 0))*Vau1)
      V0        <- (1 + sum(BqfVec))*diag(as.numeric(Vars))
    }
    Ct        <- t(ctf)%*%ctf
    
    S         <- t(Y)%*%Ct%*%Y
    CV        <- Ct%*%V0
    e         <- sum(diag(CV))
    psi       <- 2*(sum(diag(CV%*%CV)))
    k         <- psi/(2*e)
    nu        <- 2*(e^2)/psi
    Xobs      <- S/k
    pval      <- 1-pchisq(Xobs, df = nu)
    
    X2              <- as.numeric(Xobs)
    names(X2)       <- 'Chi-squared'
    nullval         <- 0
    if(test.flag[f] == 'se'){
      names(nullval)  <- 's(x)'
    } else{
      names(nullval)  <- 'g(t)'
    }
    
    comp      <- list(S = S, Ctf = Ctf, Ct = Ct, V0 = V0, ct = ctf,
                      Sigmaqf = Sigmaqf)
    out[[f]]  <- list(statistic = X2, p.value = as.numeric(pval), 
                      parameters = list(df = as.numeric(nu)), 
                      method = 'One Sample Variational ZLS Test', 
                      alternative = 'two-sided', null.value = nullval,
                      data.name = model$data.name[[f]],
                      k = as.numeric(k), psi = as.numeric(psi), e = as.numeric(e), 
                      comp = comp)
    class(out[[f]])   <- 'htest'
  }
  
  if(length(out) > 1){
    if(is.null(test)){
      return(out)
    } else{
      return(out[[test]])
    }
  } else{
    return(out[[1]])
  }
}

#### Class functions ####
plot.vsof  <- function(model, which = 1:(model$nfun + 1),
                      ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                      add.points = TRUE, xlab = NULL, ylab = NULL,
                      ...){
  oldPar  <- par()
  oldAsk  <- oldPar$ask
  par(ask = ask)
  test.flag   <- model$model$test.flag
  if(!is.null(model$model$se)){
    mf      <- model.frame(formula = model$model$fixed, data = model$model$data)
    mse     <- as.formula(paste(deparse(model$model$se),"- 1"))
    Zl      <- model.matrix(mse, data = model$model$data)
    se.dn   <- names(model.frame(formula = model$model$se, data = model$model$data))
  }
  if(!is.null(model$model$fe)){
    fe.dn  <- names(model.frame(formula = model$model$fe, data = model$model$data))
  }
  if(length(which) > 1){
    fe.counter  <- 0
    for(i in 1:model$nfun){
      if(test.flag[i] == 'se'){
        sRange    <- model$se.spline$sRange[[i]]
        if(is.null(xlab)){
          xlab  <- se.dn[i]
        }
        if(is.null(ylab)){
          ylab  <- names(mf)[1]
        }
        plot(sRange, model$model$family$linkinv(model$smoothEffect[[i]] + model$fixedEffects[1]), type = 'l', 
             ylim = range(model$model$Y), ylab = ylab, xlab = xlab, ...)
        if(add.points){
          points(Zl[,i], model$model$Y)
        }
        xlab <- NULL
      } else{
        fe.counter <- fe.counter + 1
        if(is.null(xlab)){
          xlab  <- 'time'
        }
        if(is.null(ylab)){
          ylab  <- fe.dn[fe.counter]
        }
        plot(model$model$family$linkinv(model$funcEffect[[fe.counter]]), type = 'l', xlab = xlab, ylab = ylab, ...)
        ylab <- NULL
      }
    }
    plot(model$algOut$logpd$lDelta[1:model$algOut$totalIter], xlab = 'iteration',
         ylab = '', type = 'b', pch = 16, ...)
    title(ylab = expression(paste(Delta, ' ', log(P))), line = 2, ...)
    if(ask){
      par(ask = oldAsk)
    }
  } else{
    if(which > model$nfun){
      plot(model$algOut$logpd$lDelta[1:model$algOut$totalIter], xlab = 'iteration',
           ylab = '', type = 'b', pch = 16, ...)
      title(ylab = expression(paste(Delta, ' ', log(P))), line = 2, ...)
    } else{
      if(test.flag[which] == 'se'){
        sRange    <- model$se.spline$sRange[[which]]
        if(is.null(xlab)){
          xlab  <- se.dn[which]
        }
        if(is.null(ylab)){
          ylab  <- names(mf)[1]
        }
        plot(sRange, model$model$family$linkinv(model$smoothEffect[[which]] + model$fixedEffects[1]), type = 'l', 
             ylim = range(model$model$Y), ylab = ylab, xlab = xlab, ...)
        if(add.points){
          points(Zl[,which], model$model$Y)
        }        
      } else{
        if(is.null(xlab)){
          xlab  <- 'time'
        }
        if(is.null(ylab)){
          ylab  <- fe.dn[which]
        }
        plot(model$model$family$linkinv(model$funcEffect[[which]]), type = 'l', xlab = xlab, ylab = ylab, ...)
      }
    }
  }
}

coef.vsof  <- function(model){
  model$fixedEffects
}


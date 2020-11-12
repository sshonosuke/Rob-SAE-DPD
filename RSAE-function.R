###-------------------------------------------------------------###
###    Robust SAE using density power divergence (DPD)          ###
###-------------------------------------------------------------###
library(sae)

## This script includes the following three functions: 
# 'RFH.DPD'
# 'select.alpha' 
# 'mse.RFH.DPD' 



###   Robust SAE using DPD with fixed tuning parameter 'alpha'   ###
## Input: 
# Y: m-dimensional vector of direct estimators
# X: (m,p)-matrix of covariates including an intercept
# Di: m-dimensional vector of sampling variances
# alpha: tuning parameter of DPD ('alpha=0' reduces to standard FH)
## Output:
# EB: robust empirical Bayes estimator
# beta: estimated regression coefficients
# A: estimated random effect variance 
# alpha: tuning parameter 
# itr: number of iteration 

RFH.DPD <- function(Y, X, Di, alpha=0.2, maxitr=100){
  Ep <- 10^(-5)   # torelance rate 
  
  ## weight function 
  Si <- function(Beta, A){
    mu <- as.vector(X%*%Beta)
    dnorm(Y, mu, sqrt(A+Di))^(alpha/2)
  }
  
  ## initial values
  init.fit <- eblupFH(Y~X[,-1], vardir=Di)
  hBeta <- init.fit$fit$estcoef$beta
  hA <- init.fit$fit$refvar
  if(hA<0.1){ hA <- 0.1 }
  
  ## iteration method 
  for(k in 1:maxitr){
    c.hBeta <- hBeta 
    c.hA <- hA
    # update Beta
    W <- diag(Si(hBeta, hA)/(hA+Di))
    hBeta <- as.vector( solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%Y )
    # update A
    ss <- Si(hBeta, hA)
    hmu <- as.vector(X%*%hBeta)
    V <- 1/sqrt(2*pi*(hA+Di))
    vv <- ss - alpha/(alpha+1)^(3/2) * V^alpha
    num <- ( (Y-hmu)^2*ss - Di*vv ) / (hA+Di)^2
    denom <- vv / (hA+Di)^2
    hA <- max(0.001, mean(num)/mean(denom))
    # convergence check 
    dd <- 100*sum(abs(hBeta - c.hBeta))/sum(abs(c.hBeta))
    if( dd<Ep ){ break }
  }
  
  ## EB estimator
  hmu <- as.vector(X%*%hBeta)
  sh <- Si(hBeta, hA)
  EB <- Y - sh*Di*(Y-hmu)/(hA+Di)
  
  Res <- list(EB=EB, beta=hBeta, A=hA, alpha=alpha, itr=k)
  return(Res)
}


###   tuning parameter selection via excess MSE    ###
## Input: 
# Y: m-dimensional vector of direct estimators
# X: (m,p)-matrix of covariates including an intercept
# Di: m-dimensional vector of sampling variances
# IR: percentage of inflation rate of MSE (default: 5%)
## Output:
# selected alpha 

select.alpha <- function(Y, X, Di, IR=5, maxitr=100){
  ## function for excess MSE
  EMSE <- function(a){
    fit <- RFH.DPD(Y, X, Di, alpha=a)
    mu <- as.vector( X%*%fit$beta )
    A <- fit$A
    V <- 1/sqrt(2*pi*(A+Di))
    vv <- V^(2*a)/(2*a+1)^(3/2) - 2*V^a/(a+1)^(3/2) + 1
    rate <- 100* mean(vv*Di^2/(A+Di)) / mean(A*Di/(A+Di)) - IR
    return(rate)
  }
  
  ## bi-section method 
  Ep <- 10^(-4)
  lw <- 0.0001
  up <- 0.5
  alpha <- (lw+up)/2
  v1 <- EMSE(lw)
  v2 <- EMSE(up)
  for(k in 1:maxitr){
    c.alpha <- alpha
    vc <- EMSE(alpha)
    if( vc*v1<0 ){ 
      up <- alpha
      alpha <- (lw+up)/2
      v2 <- vc
    }
    if( vc*v2<0 ){
      lw <- alpha
      alpha <- (lw+up)/2
      v1 <- vc
    }
    if( abs(alpha-c.alpha)<Ep ){ break }
  }
  
  ## output 
  return(alpha)
}





###   MSE estimation via parametric boostrap    ###
## Input: 
# Y: m-dimensional vector of direct estimators
# X: (m,p)-matrix of covariates including an intercept
# Di: m-dimensional vector of sampling variances
# alpha: tuning parameter of DPD ('alpha=0' reduces to standard FH)
# B: number of bootstrap samples (default: 200)
## Output:
# MSE: second-order unbiased MSE 
# nMSE: first-order unbiased MSE

mse.RFH.DPD <- function(Y, X, Di, alpha=0.2, B=200){
  m <- length(Y)
  
  ## weight function 
  Si <- function(Beta, A){
    mu <- as.vector(X%*%Beta)
    dnorm(Y, mu, sqrt(A+Di))^(alpha/2)
  }
  
  ## function for MSE (first order)
  g12 <- function(Beta, A){
    mu <- as.vector(X%*%Beta)
    V <- 1/sqrt(2*pi*(A+Di))
    vv <- (2*alpha+1)^(-3/2)*V^(2*alpha) - 2*(alpha+1)^(-3/2)*V^(alpha) + 1
    return( A*Di/(A+Di) + Di^2/(A+Di)*vv )
  }
  
  ## fit 
  fit <- RFH.DPD(Y, X, Di, alpha=alpha)
  Beta <- fit$beta
  A <- fit$A
  sh <- Si(Beta, A)
  mu <- as.vector( X%*%Beta )
  V <- 1/sqrt(2*pi*(A+Di))
  
  ## quantities for asymptotic variance
  Jb <- t(X)%*%diag(V^alpha/(A+Di))%*%X/(alpha+1)^(3/2)/m
  Ja <- mean( V^alpha/(A+Di)^2*(2-alpha)*(alpha^2+alpha+1)/(alpha+1)^(5/2)/2 )
  Kb <- t(X)%*%diag(V^(2*alpha)/(A+Di))%*%X/(2*alpha+1)^(3/2)/m
  h <- 2*(2*alpha^2+1)/(2*alpha+1)^(5/2)-alpha^2/(alpha+1)^3
  Ka <- mean( V^(2*alpha)/(A+Di)^2*h )
  
  ## quantities for MSE
  g3 <- Di^2*V^(2*alpha)/(A+Di)^2/(2*alpha+1)^(3/2)*X%*%solve(Jb)%*%Kb%*%solve(Jb)%*%t(X)
  g3 <- diag(g3)
  g4 <- Di^2*V^(2*alpha)/(A+Di)^3/(2*alpha+1)^(7/2)*(alpha^4-0.5*alpha^2+1)*Ka/Ja^2
  
  ## parametric bootstrap 
  boot.g12 <- matrix(NA, B, m)
  boot.g5 <- matrix(NA, B, m)
  for(b in 1:B){
    boot.Y <- mu + rnorm(m, 0, sqrt(A+Di))
    boot.fit <- RFH.DPD(Y=boot.Y, X=X, Di=Di, alpha=alpha)
    boot.Beta <- boot.fit$beta
    boot.A <- boot.fit$A
    boot.g12[b,] <- g12(boot.Beta, boot.A)   # first order terms
    boot.EB <- boot.fit$EB
    boot.s <- Si(boot.Beta, boot.A)
    EB1 <- boot.Y - Di/(A+Di)*sh*(boot.Y-mu)
    EB2 <- boot.Y - Di/(A+Di)*(boot.Y-mu)
    boot.g5[b,] <- 2*(boot.EB-EB1)*(EB1-EB2)
  }
  
  ## summary 
  mse1 <- g12(Beta, A)
  mse2 <- 2*mse1 - apply(boot.g12, 2, mean) + g3/m + g4/m + apply(boot.g5, 2, mean)
  Res <- list(MSE=mse2, nMSE=mse1)
  return(Res)
}




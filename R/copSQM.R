
# This is an  function named 'Copula based Stochastic frontier Quantile model'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#

ptALD <- function(u, sigu,tau){
  ptrunc(u, "alaplace", a=0, b=Inf,location=0, scale=sigu, kappa=tau)
}

rtALD <- function(n, mean, sigu,tau){
  rtrunc(n, "alaplace", a=0, b=Inf,location=mean, scale=sigu, kappa=tau)
}



## Main Function
copSQM=function(Y,X,family,tau,RHO,LB,UB){

  like<-function(theta,Y,X,family,tau,RHO,LB,UB){


    midX=X
    midY=Y

    XX=as.matrix(cbind(1,midX))
    K=ncol(XX)
    sigmav=abs(theta[K+1])
    sigmau=abs(theta[K+2])
    rho=theta[K+3]
    n=length(midY)
    m=n


    w=c(midY-t(theta[1:K])%*%t(XX))
    set.seed(1988)


    u=abs(replicate(n,rALD(n,0,sigmau,tau)))
    W=t(replicate(n,w))
    gv =  dALD(u+W,mu=0,sigma=sigmav,p=tau)+0.000001

    gv=matrix(t(gv),nrow=n,ncol=n)
    Gv=pALD(u+W,mu=0,sigma=sigmav,p = tau )
    Fu=ptALD(u,sigmau,tau)
    Gvv=c(Gv)
    mm=length(Fu)
    for ( i in 1:mm){
      if (is.infinite(Fu[i]))  # control for optimization
        Fu=0.0000000000000001
      if (is.infinite(Gvv[i]))  # control for optimization
        Gvv=0.000000000000001
      if (is.nan(Fu[i]))  # control for optimization
        Fu=0.00000000000001
      if (is.nan(Gvv[i]))  # control for optimization
        Gvv=0.00000000000001
    }

    if (family==2){
      rho2=theta[length(theta)]
      gaucopula=BiCopPDF(Fu, Gvv, family=family, par=rho, par2=rho2)+0.00000001
    }else{
      gaucopula=BiCopPDF(Fu, Gvv, family=family, par=rho, par2=0)+0.00000001}

    gaucopula=matrix(gaucopula,nrow=m,ncol=n)
    hw=sum(log(1/m*diag(gv%*%gaucopula)))
    if (is.infinite(hw))  # control for optimization
      hw<--n*100
    if (is.nan(hw))  # control for optimization
      hw<--n*100

    cat("Sum of log Likelihood for CSQM ->",sprintf("%4.4f",hw),"\n")


    return(hw) # log likelihood

  }

  ### End Function #############3
  #=================================================


  # start here
  cc=coef(lm(Y~X))
  K=length(cc)
  if (family==2){
    lower =c(rep(-Inf,K),0.01,0.01,LB,4.1)
    upper =c(rep(Inf,K+2),UB,50)

    theta=c(cc,sigmav=1,sigmau=1,rho=RHO,df=4)
  }else{
    lower =c(rep(-Inf,K),0.01,0.01,LB)
    upper =c(rep(Inf,K+2),UB)
    theta=c(cc,sigmav=1,sigmau=1,rho=RHO)
  }


  model <- optim(theta,like,Y=Y,X=X,family=family,tau=tau,
                 control = list(maxit=100000,fnscale=-1),method="L-BFGS-B",
                 lower =lower,upper =upper, hessian=TRUE )
  # table of results
  coef<- model$par
  k=length(coef)
  model$se <- sqrt(-diag(solve(model$hessian)))

  for(i in 1:k){
    if (is.nan(model$se[i]))  # control for optimization
      model$se[i] <- sqrt(-diag(solve(-model$hessian)))[i]
  }

  n=length(Y)
  S.E.= model$se
  (paramsWithTs = cbind (model$par , coef/S.E. ) )
  stat=coef/S.E.
  pvalue <- 2*(1 - pnorm(abs(stat)))
  result <- cbind(coef,S.E.,stat,pvalue)
  result
  BIC= -2*model$value+ (log(n)*length(coef))
  AIC = -2*model$value + 2*length(coef)


  output=list(
    result=result,
    AIC=AIC,
    BIC=BIC,
    Loglikelihood=model$value
  )
  output

}


## Required packages
#library(truncnorm)
#library(mvtnorm)
#library("VineCopula")
#library("frontier")
#library(ald)
#library("LaplacesDemon")

# example included in FRONTIER 4.1 (cross-section data)
#data(front41Data)
#attach(front41Data)
# Cobb-Douglas production frontier
#cobbDouglas <- sfa( log(output)~log(capital)+log(labour),data=front41Data)
#summary(cobbDouglas)


# Select familty  copula upper and lower bouubd ( look at Vinecopula package)
#family=1   # 1 is Gaussian, 2 is Student-t, 3 is Clayton and so on....

#Gaussian (-.99, .99)
#Student t (-.99, .99)
#Clayton (0.1, Inf)
#Y=log(output)
#X=cbind(log(capital),log(labour))
#model=copSQM(Y=Y,X=X,family=1,tau=0.5,RHO=0.5,LB=-0.99,UB=0.99)


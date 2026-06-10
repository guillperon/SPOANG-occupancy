## Code appendix to "Residuals of species distribution model reveal hotspots and strongholds for a passerine affected by the cagebird trade "
## Author : Guillaume Péron

## Negative Log-likelihood function    /!\ uses specific covariate names  "var1" to "var7"
## -----------------------------------------------------------------------------
Loglik = function (theta) {
  # parameters for initial occupancy
  beta.psi0 = theta[1:npsi0] ; i=npsi0
  # parameters for detection
  beta.P = theta[i+1:nP] ; i=i+nP
  if(doPtime) {
  beta.Pt = theta[i+1] ; i=i+1
  beta.Pt2 = theta[i+1] ; i=i+1
  } else { beta.Pt=beta.Pt2=0 }
  if(doPinter) {
  beta.PtxZ = theta[i+1] ; i=i+1
  } else { beta.PtxZ=0 }
  # parameters for extinction
  beta.ext = theta[i+(1:nex)] ; i=i+nex  
  # parameters for colonization
  beta.col = theta[i+(1:ncol)] ; i=i+ncol  
  # parameters for spatial autoregression effects (with transformation to force direction of effects as implemented in the main text) 
  if(doALext) { beta.AL.ext = -3*plogis(theta[i+1]) ; i=i+1  }  else  beta.AL.ext = 0
  if(doALcol) { beta.AL.col = +3*plogis(theta[i+1]) ; i=i+1  }  else  beta.AL.col = 0
  
  names(beta.psi0) = covarid.psi0
  names(beta.P) = covarid.P
  names(beta.ext) = covarid.ext
  names(beta.col) = covarid.col
     
  ## Year 1: initial state
  PSIupdate = PSI = plogis( COVAR[,covarid.psi0] %*% (beta.psi0) )
  P = plogis( COVAR[,covarid.P] %*% (beta.P) + beta.Pt*tt[1] + beta.Pt2*tt2[1] + beta.PtxZ*tt[1]*COVAR[,"var7"] )
  LL=0
  for ( k in 1:Ncells ) {
    psi=PSI[k]
    p = P[k]
    nn = nvisits[k,1]
    if(nn>0) {
    mm = nobs[k,1]
    likk = (mm==0)*(1-psi + psi*(1-p)^nn) + (mm>0)*(psi*dbinom(mm,nn,p))
    LL=LL - log(likk)
    }
  }

  # Years 2 to T : extinctions and colonisations
  for (year in 2:Nyears) {
    P = plogis( COVAR[,covarid.P] %*% (beta.P) + beta.Pt*tt[year] + beta.Pt2*tt2[year] + beta.PtxZ*tt[year]*COVAR[,"var7"] )
    for ( k in 1:Ncells ) {
      p = P[k]
      psi = PSI[k]
      # Average occupancy in the neighborhood  
      psi_bar= mean(PSI[NEIGH[[k]]],na.rm=TRUE) 
      if(is.na(psi_bar)) psi_bar=0.001  # no neighbors in the dataset. In practice, these were small islands. Otherwise set to mean(PSI) or smth else
      # Extinction and colonization probabilities
      ext = plogis( sum(beta.ext*COVAR[k,covarid.ext]) + beta.AL.ext*psi_bar )
      col = plogis( sum(beta.col*COVAR[k,covarid.col]) + beta.AL.col*psi_bar )
      # Update okupancy
      psiup =  psi*(1-ext) + (1-psi)*col  
      PSIupdate[k] = psiup
      # Likelihood of data bit
      nn = nvisits[k,year]
      if(nn>0) {
      mm = nobs[k,year]
      likk = (mm==0)*(1-psi + psi*(1-p)^nn) + (mm>0)*(psi*dbinom(mm,nn,p))
      LL=LL - log(likk)
      }
    }
    PSI <- PSIupdate
  }
    return(LL)
}
    
## Test on a simplistic scenario in which most estimates should be zero
## -----------------------------------------------------------------------------
## Simulate dummy occupancy data
simulation0 <- function(psi, p, n, K) {
T = length(psi)
kmax = max(K)
ans = matrix(NA, n, T)
if (T != length(p)) cat("Wrong entry\n")
if (T != dim(K)[2]) cat("Wrong entry\n")
if (n != dim(K)[1]) cat("Wrong entry\n")
for (i in 1:n) {
     for (t in 1:T) {
          ans[i,t] = rbinom(1,1,psi[t])*rbinom(1,K[i,t],p[t])
     }
}
return(ans)
}
Ncells = 100  ## number of cells
Nyears=10     ## number of years
K =10         ## number of visits per year and cell
nvisits = matrix(K,Ncells,Nyears)
nobs = simulation0(psi=seq(0.9,0.6,length.out=Nyears),p=seq(0.6,0.9,length.out=Nyears),n=Ncells,K=nvisits)

## Simulate dummy covariate data
COVAR=matrix(rnorm(Ncells*7,0,1),Ncells,7)
colnames(COVAR) = paste("var",1:7,sep="")

## Create random reighborhoods
NEIGH = vector(mode = "list", length = Ncells)
for(i in 1:Ncells) {
    NEIGH[[i]] = unique(c(NEIGH[[i]], sample(1:Ncells, max(0, 8-length(NEIGH[[i]])), replace=F)))
    for(j in NEIGH[[i]]) NEIGH[[j]]=unique(c(NEIGH[[j]],i))
}

## Cammands to specify the model (same structure as "Full model" in the main text)
 covarid.psi0=colnames(COVAR)[1:7]  ## names of the covariates to include in the component for initial occupancy
 covarid.P=colnames(COVAR)[1:7]     ## names of the covariates to include in the component for detection
 covarid.ext=colnames(COVAR)[1:7]   ## names of the covariates to include in the component for extinction
 covarid.col=colnames(COVAR)[1:7]   ## names of the covariates to include in the component for colonization
 doALext=TRUE                      ## do neighborhood rescue effect on extinction
 doALcol=TRUE                      ## do neighborhood rescue effect on colonisation
 doPtime=TRUE                      ## do time variation in P
 doPinter=TRUE                     ## do interaction between time and covariate #7 on P
 npsi0 = length(covarid.psi0)
 nP = length(covarid.P)
 nex = length(covarid.ext)
 ncol = length(covarid.col)
tt = ((1:Nyears)-mean(1:Nyears))/sd(1:Nyears)  ## time variable
tt2 = tt^2    

## Run optimization
theta_init = c(qlogis(0.9),rep(0,6), qlogis(0.6),rep(0,6), 0.1,0,0, qlogis(0.1),rep(0,6), qlogis(0.01),rep(0,6), 0,0)
ss <- nlminb(theta_init, Loglik, control=list(iter.max=100,rel.tol=1e-9, x.tol=1e-7), lower = -5, upper = 5)    

## Compute standard errors 
## /!\ Might generate an error msg 
## which means either optimization was unsuccessful or there are boundary parameters 
## Try drop/fix boundary parameters, or different initial values
library(numDeriv)
hess <- hessian(Loglik, ss$par)
std = sqrt(diag(solve(hess)))    

## A few checks
plot(ss$par~theta_init)   ## should have 5 points on diagonal, all the rest is noise
c(ss$par[15], std[15])    ## should be >0 if model catches the increase in detection
c(plogis(ss$par[18]), plogis(ss$par[18])/(1-plogis(ss$par[18]))*std[18])    ## should be >0 if model catches the decrease in occupancy
c(plogis(ss$par[25]), plogis(ss$par[25])/(1-plogis(ss$par[25]))*std[25])    ## should be near zero

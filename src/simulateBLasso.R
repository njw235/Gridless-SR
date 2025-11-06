simulateBLasso = function(M, seed0, warmup=10000){
  
    # load data
    load("SupportReduction/julia/concreteStrength.Rdata")
    # load bayesLasso sampler
    source("SupportReduction/julia/bayesLasso.R")
    # load "ground truth" 
    load("SupportReduction/julia/concrete_truth.Rdata")
    
    # sample
    set.seed(seed0)
    ch = run.bl(X=X_concrete_scaled,Y=Y_concrete-mean(Y_concrete),lambda=50,fast=TRUE,M=M,burnin=warmup)
    ch = cbind(ch$beta,ch$sigma2)
    colnames(ch)=c(colnames(X_concrete_scaled),"sigma2")
    
    meanTruth =concrete_truth$meanTruth
    varTruth = concrete_truth$varTruth
    phiTruth = concrete_truth$phiTruth
    
  
  
  out = list(x=ch,varTruth=varTruth,meanTruth=meanTruth,phiTruth=phiTruth)
  return(out)
}

### Example ###
# simulate length M chain for each beta_j, sigma
ch = simulateBLasso(M=1000,seed0 = 12345,warmup = 10000)

# check
library(momentLS)
diag(ch$varTruth)
sapply(1:ncol(ch$x), function(j) avar(SR1(autocov(ch$x[,j]),delta = tune_delta(ch$x[,j],nSplits = 5)$delta*0.8)))

library(momentLS)

K = 50 # number of discrete states

# simulate g, pi, Q (should be fixed for all simulations to use the same dynamic)
set.seed(1356)
g = matrix(rnorm(K*1),ncol=1)
discreteMC = simulate_discreteMC(nStates = K)

# Simulate a finite state space MH chain
set.seed(1359) # vary this seed to get different chains
chainParams = list(type="MH",  M = 10000,  nStates = K, g = g, discreteMC = discreteMC, d = 1)
ch_mh = generateChain(chainParams)
r = autocov(ch_mh$x)

# check
ch_mh$varTruth
avar(SR1(r,delta = tune_delta(ch_mh$x,nSplits = 5)$delta*0.8))

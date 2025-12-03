include("SR1.jl")
R"library(momentLS)"

R"source('simulateMH.R')"
errors = zeros(51)
terrors = zeros(51)
R"measlist = vector('list', length = 50)"
R"set.seed(1234)"
R"g = matrix(rnorm(50), ncol = 1)"
R"discreteMC = simulate_discreteMC(nStates=50)"
R"chainParams = list(type = 'MH', M = 1000, nStates = 50, g = g, discreteMC = discreteMC, d = 1)"
for j in 1:51

    R"ch = generateChain(chainParams)"

    R"r = autocov(ch$x)"
    R"dhat = tune_delta(ch$x,5)$delta*0.8"
    @rget r
    @rget dhat
    m = SR1_gridless(r, dhat,"unweighted", [0.0], [0.0], 1e-9)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight

    R"err = L2diff_L2Moment(r, supp, weight)"
    R"terr = L2diff_twoMoments(supp, weight, ch$F_support, ch$F_weights)"
    @rget terr
    @rget err
    errors[j] = err
    terrors[j] = terr
    R"measlist[[j]] = cbind(supp, weight)"
end
R"saveRDS(measlist, 'MHmeas.RDS')"
println(terrors)

open("MH.txt", "w") do io
    println(io, errors)
    println(io, terrors)
end

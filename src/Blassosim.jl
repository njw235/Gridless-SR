include("momentLS.jl")
R"load('MC_chains.Rdata')"
R"library(momentLS)"
R"print(getwd())"
R"source('simulateBLasso.R')"
R"measlist = vector('list', length = 50)"
i = ARGS[1]
errors = zeros(50)
for j in 1:50
    @rput i
    @rput j
    R"ch_blasso = simulateBLasso(M =1000, seed0 = 1234 + j, warmup = 10000)"
    R"r = autocov(ch_blasso$x[,as.integer(i)])"
    R"dhat = tune_delta(ch_blasso$x[,as.integer(i)],5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0], [0.0], 1e-9)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight

    R"err = L2diff_L2Moment(r, supp, weight)"
    @rget err
    errors[j] = err
    println(err)
    R"measlist[[j]] = cbind(supp, weight)"
end
R"print(measlist)"
R"saveRDS(measlist, paste(i, 'Blasso.rds', sep = ''))"
open(string("BLasso", i, ".txt"), "w") do io
    println(io, errors)
end

include("momentLS.jl")
R"load('MC_chains.Rdata')"
R"library(momentLS)"

errorbP = zeros(8,2)
errormhP = zeros(5,2)

for i in 1:8
    @rput i
    R"r = autocov(ch_blasso$x[,i])"
    R"dhat = tune_delta(ch_blasso$x[,i],5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0], [0.0], 3e-16)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight
    R"m = SR1(r, dhat)"

    R"err2 = L2diff_L2Moment(r, supp, weight)"
    R"err3 = (asympVariance(weight, supp) - ch_blasso$varTruth[i,i])^2"


    @rget err2
    @rget err3


    errorbP[i,1] = err2
    errorbP[i,2] = err3

end

#= for i in 1:5
    @rput i
    R"r = autocov(ch_mh$x[,i])"
    R"dhat = tune_delta(ch_mh$x[,i],5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0], [0.0], 3e-16)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight

    R"err2 = L2diff_L2Moment(r, supp, weight)"
    R"err3 = (asympVariance(weight, supp) - ch_mh$varTruth[i,i])^2"

    @rget err2
    @rget err3

    errormhP[i,1] = err2
    errormhP[i,2] = err3
end
=#

#println(errormhP)
println(errorbP)
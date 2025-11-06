include("momentLS.jl")
R"library(momentLS)"

R"source('SupportReduction/julia/simulateMH.R')"
errors = zeros(50)
measures = Matrix{Float64}
R"set.seed(1234)"
R"g = matrix(rnorm(50), ncol = 1)"
R"discreteMC = simulate_discreteMC(nStates=50)"
R"chainParams = list(type = 'MH', M = 1000, nStates = 50, g = g, discreteMC = discreteMC, d = 1)"
for j in 1:50

    R"ch = generateChain(chainParams)"

    R"r = autocov(ch$x)"
    R"dhat = tune_delta(ch$x,5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0], [0.0], 1e-9)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight
    R"m = SR1(r, dhat)"

    R"err = L2diff_L2Moment(r, m$support, m$weights) - L2diff_L2Moment(r, supp, weight)"
    @rget err
    errors[j] = err
    measure = hcat(supp, weight)
	println(measure)
	global measures = vcat(measures, measure)
end


open("MH.txt", "w") do io
    println(io, errors)
    println(io, measures)
end
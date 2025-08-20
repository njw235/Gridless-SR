include("momentLS.jl")
R"library(momentLS)"

R"set.seed(1234)"

error = zeros(3)
errora = zeros(3)
n = [50,100]
for i in 1:2
    N = n[i]
    @rput N
    R"x = generateChain(list(type  = 'AR', rho = 0.9, M = N))$x"
    R"r = autocov(x)"
    R"dhat = tune_delta(x,5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0],[0.0], 3e-8)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight
    m2 = momentLS(-1+dhat, 1-dhat, r, 3e-8)
    s2 = m2[1]
    w2 = m2[2]
    @rput s2
    @rput w2
    R"err1 = L2diff_L2Moment(r, supp, weight)"
    R"err2 = L2diff_L2Moment(r, s2, w2)"
    @rget err1
    @rget err2
    error[i] = err1
    errora[i] = err2
end

println(error)
println(errora)

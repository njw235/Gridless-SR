include("momentLS.jl")
R"library(momentLS)"

errorps = zeros(50)
error1 = zeros(50)
error2 = zeros(50)
error3 = zeros(50)
error4 = zeros(50)


R"set.seed(9999)"
R"rho = c(runif(1, -0.9, -0.7), runif(1, -0.6, -0.4), runif(1, -0.2, 0.2), runif(1, 0.4, 0.6), runif(1, 0.7, 0.9))"
@rput i
R"rh = rho[as.integer(i)]"
N = 1000
for j in 1:10
    @rput N
    R"x = generateChain(list(type  = 'AR', rho = rh, M = N))$x"
	R"r = autocov(x)"
    R"dhat = tune_delta(x,5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0],[0.0], 3e-8)
    supp = m[1]
	weight = m[2]
    @rput supp
	@rput weight
	R"m = SR1(r,dhat, gridControl = list(cm = FALSE, scale = 'equidist'))"
    R"m2 = SR1(r, dhat, n_alphas = 101,gridControl = list(cm = FALSE, scale = 'equidist'))"
	R"m3 = SR1(r, dhat, n_alphas = 51,gridControl = list(cm = FALSE, scale = 'equidist'))"
	R"m4 = SR1(r, dhat, n_alphas = 11, gridControl = list(cm= FALSE, scale = 'equidist'))"
    R"e1 = L2diff_L2Moment(r,m$support, m$weights)"
	R"e2 = L2diff_L2Moment(r, m2$support, m2$weights)"
    R"e3 = L2diff_L2Moment(r, m3$support, m3$weights)"
	R"e4 = L2diff_L2Moment(r, m4$support, m4$weights)"
	@rget e1
    @rget e2
	@rget e3
	@rget e4
    error1[j] = e1
	error2[j] = e2
    error3[j] = e3
	error4[j] = e4
    R"errorp = L2diff_L2Moment(r,supp, weight)"
	@rget errorp
	errorps[j] = errorp
end

@rget rh


ep1 = error1 - errorps
ep2 = error2 - errorps
ep3 = error3 - errorps
ep4 = error4 - errorps


d1 = mean(ep1)
d2 = mean(ep2)
d3 = mean(ep3)
d4 = mean(ep4)

v1 = var(ep1)
v2 = var(ep2)
v3 = var(ep3)
v4 = var(ep4)

open(string("MCMCout",ARGS[1],".txt"), "w") do io
    println(io, d1)
    println(io, d2)
    println(io,d3)
	println(io, d4)
    println(io, v1)
    println(io, v2)
    println(io,v3)
	println(io, v4)
end
#N = 6 .^[2:1:5;]
#plot(N, [d1 d2 d3], yerr = [v1 v2 v3], label = [501 101 51], xaxis=:log)
#title!(L"$\rho = %$rh$")
#ylabel!(L"Difference in $\ell_2$ error")
#xlabel!(L"Sample Size")

#savefig(string("plot", ARGS[1], ".png"))



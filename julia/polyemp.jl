#include("momentLS.jl")

using DynamicPolynomials
using SumOfSquares
using Random
using Clarabel
using Distributions

transform = function(x,i, J, N)
    if(i < J)
	    return(sign(i)*((1+x)*(1-2.0^-abs(i)) - x))
    else
        return(1 - 1/N - x/N)
    end
end

estimate_poly = function(i,r)
	m = Int(ceil(exp(1+1/exp(1))*log(10^12)))
	t = Int(floor(2^abs(i) * log(10^12)))
	if i > ceil(log(length(r)))
		a0 = 1-1/length(r)
	else
		a0 = (1- 2.0^-abs(i))
	end
	up = min(m-1,t)

	b = zeros(up)

	for j in 0:up-1
		for k in j:min(t,length(r)-1)
			b[j+1] += sign(i)^k * r[k+1] * binomial(big(k), big(j)) * a0^(k-j)* (1-a0)^j
		end
		b[j+1] = (-1)^j * b[j+1]
	end
	return(b,a0,up)
end


solve_LDA = function(f,g)
    
    J = Int(ceil(log(length(f))))
    optim = zeros(J)
    id = [1:1:J;]
	p = Dict(id .=> [estimate_poly(i,f) for i in [1:1:J;]])
    time = 0
    i = 1
    for ind in id
        model = SOSModel(Clarabel.Optimizer)
        set_string_names_on_creation(model, false)
		@polyvar x
		fx = sum(p[ind][1] .* transform(x, ind, J, length(f)) .^[0:1:(length(p[ind][1])-1);])
        gx = sum(g .* transform(x,ind,J, length(f)) .^ [0:1:4;])

        S = @set x >= 0 && 1-x >= 0

        @variable(model,s)
			
		@constraint(model,c, fx >= s*gx, domain = S)
		@objective(model, Max, s)
		optimize!(model)
        time += solve_time(model)
        optim[i] = objective_value(model)
    end
    min = findmin(optim)[1]
    return(time, min)
end

times = zeros(length(10:25:250;))
LDAtimes = zeros(length(10:25:250;))
val = zeros(length(10:25:250;))
LDAval = zeros(length(10:25:250;))
for i in [10:25:250;]

    @polyvar x
    f_c = rand(Uniform(0,1), i+1)
    g_c = rand(Uniform(0,1), 5)
    f = sum(f_c .* x.^[0:1:i;])
    g = sum(g_c .* x.^[0:1:4;])
    S = @set x>=0 && x <= 1
    model = SOSModel(Clarabel.Optimizer)
    @variable(model, s)
    @objective(model, Max, s)
    @constraint(model, c, f -s*g>= 0,domain = S)
    optimize!(model)
    id = Int((i-10)/25) + 1
    times[id] = solve_time(model)
    val[id] = objective_value(model)
    LDAtimes[id], LDAval[id] = solve_LDA(f_c,g_c)
end

println(times)
#plot([10:25:250;], [times LDAtimes])
#savefig("timeplt.png")
#plot([10:25:250;], abs(val-LDAval))

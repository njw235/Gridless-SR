using Clarabel
using SumOfSquares
using DynamicPolynomials
using Random
using Distributions
using LinearAlgebra
using LinearSolve
using Plots
using StatsBase


transform = function(x,i)
	return(sign(i)*((1+x)*(1-2.0^-abs(i)) - x))
end


Empirical = function(x)
    
    n = length(x)
    supp = [0:1:findmax(x)[1];]
    count = proportions(x)
    
    return(count)

end

grad_opt = function(p, supp, weight, solver,delta,x)
	gradients = zeros(13)
	supports = zeros(13)
	n = -Int(floor(log2(delta)))
	for ind in zip([1:1:n;],[1:1:n;])
		#trying with replacing alpha with x
		model = SOSModel(solver)
        set_string_names_on_creation(model, false)
		f = 0
		for i in 1:length(supp)
			g=1
			for j in 1:length(supp)
                if(j != i)
					g = g* (1 - transform(x,ind[2])*supp[j])
                end
			end
			f = f + weight[i]*g
		end
		d = prod((1 .- transform(x,ind[2]).*supp))
	
			
		h = (f - p[ind[2]]*d)*(1-transform(x,ind[2]))

		S = @set x >= 0 && 1-x >= 0

		@variable(model,s)
			
		@constraint(model,c, h >= s*d, domain = S)
		@objective(model, Max, s)
		set_silent(model)
		optimize!(model)
			
		v = moment_matrix(model[:c])
		pt = atomic_measure(v, 1e-4)
		if(typeof(pt) != Nothing)
			if(length(pt.atoms[1].center) == 1)
				supports[ind[1]] = transform(pt.atoms[1].center[1],ind[2])
				
			else
				supports[ind[1]] = transform(pt.atoms[1].center[2], ind[2])
			end
		end
		if(is_solved_and_feasible(model))
			gradients[ind[1]] = objective_value(model)
		else
			gradients[ind[1]] = 1000
		end
	end
	return(gradients, supports)
end

SRm = function(supp, weight,r)
    validmeasure = false
    support = copy(supp)
    current = copy(weight)
    exponent = [0:1:length(r)-1;]
    while(!validmeasure)
        B = zeros(length(support), length(support))
        for i in 1:length(support)
            for j in 1:length(support)
                B[i,j] = 1/(1-support[i]*support[j])	
            end
        end
        c = zeros(length(support))
        for i in 1:length(support)
            c[i] = sum((support[i].^exponent) .* r)
        end

        prob = LinearProblem(B,c)
        sol = LinearSolve.solve(prob)
        new = sol.u
        
        if(all( >=(0), new))
            validmeasure = true
            weight = new
        else
            t = - current ./(new .- current)
            pop!(t)
            t[t .< 0] .= typemax(Int)
            t[t .> 1] .= typemax(Int)
            bd = findmin(t)
            current = (1-bd[1]).*current + bd[1] .* new
            deleteat!(support, bd[2])
            deleteat!(current,bd[2])
        end

    end
return(support, weight)
end

estimate_poly = function(i,r,x)
	m = Int(ceil(exp(1+1/exp(1))*log(10^6)))
	t = Int(floor(2^abs(i) * log(10^6)))
	a0 = (1- 2.0^-abs(i))
	up = min(m-1,t)

	b = zeros(up)

	for j in 0:up-1
		for k in j:min(t,length(r)-1)
			b[j+1] += sign(i)^k * r[k+1] * binomial(big(k), big(j)) * a0^(k-j)* (1-a0)^j
		end
		b[j+1] = (-1)^j * b[j+1]
	end
	
    return(sum(b .* x .^ [0:1:up-1;]))
end

mixingmeasure = function(r, delta,supp, weight, tol, graph = false)
    @polyvar x
	id = [1:1:13;]
	dictionary = Dict(id .=> [estimate_poly(i,r,x) for i in [1:1:13;]])
	pts = [0:0.01:1-delta;]
	solver = Clarabel.Optimizer
	n = length(r)
	exponents = [0:1:n-1;]
    s = copy(supp)
    w = copy(weight)
	conv = false
	count = 0
	while(count < 200 && !conv)
		SRstep = SRm(s, w,r)
		s = SRstep[1]
		w = SRstep[2]
	
		
		points = grad_opt(dictionary, s, w, solver,delta,x)
		index = findmin(points[1])[2]
		if(findmin(points[1])[1] > -tol)
			conv = true
		end
		append!(s, points[2][index])
		append!(w, 0)
		if(graph == true)
            a = zeros(length(pts))
            b = zeros(length(pts))
            for i in 1:length(a)
                a[i] = -sum(r.* pts[i].^exponents)
                b[i] = sum(weight.*(1 .- supp)./(1 .- pts[i].*supp))
            end
            val = (1 .- pts) .* (a+b)
            if(count == 1)
                display(plot(pts, val))
            else
                display(plot!(pts,val))
            end
        end
            
		count = count + 1
	end
	if(!conv)
        print("failed to converge")
    end
	return(s, w)
end




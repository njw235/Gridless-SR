import Clarabel
using SumOfSquares
using DynamicPolynomials
using RCall
using Random
using Distributions
using LinearAlgebra
using LinearSolve
using Plots
using HomotopyContinuation


transform = function(x,i)
	return(sign(i)*((1+x)*(1-2.0^-abs(i)) - x))
end


grad_optimize = function(r,p, supp, weight,delta)
	
	solver = Clarabel.Optimizer
	gradients = zeros(26)
	supports = zeros(26)
	n = -Int(floor(log2(delta)))
	for ind in zip([1:1:2*n;],append!([1:1:n;], [-n:1:-1;]))
		#trying with replacing alpha with x
		model = SOSModel(solver)
        set_string_names_on_creation(model, false)
		@polyvar x
		f = 0
		for i in 1:length(supp)
			g=1
			for j in 1:length(supp)
				if( j != i)
					g = g* (1 - transform(x,ind[2])*supp[j])
				end
			end
			f = f + weight[i]*(1 + transform(x, ind[2])*supp[i])*g
		end
		d = prod((1 .- transform(x,ind[2]).*supp))
	
			
		h = f + (-2 * sum(p[ind[2]][1] .* x.^[0:1:p[ind[2]][3]-1;]) + r[1])*d
		if abs(ind[2]) == n
			S = @set x>= 2^(abs(n)) * delta - 1 && 1-x >= 0
		else
			S = @set x >= 0 && 1-x >= 0
		end
		@variable(model,s)
			
		@constraint(model,c, h >= s*d, domain = S)
		@objective(model, Max, s)
		set_silent(model)
		optimize!(model)
		sol = SemialgebraicSetsHCSolver(; compile = false)	
		v = moment_matrix(model[:c])
		pt = atomic_measure(v, 1e-4, sol)
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
momentLS = function(a, b, r, tol)
	pts = []
	solver = Clarabel.Optimizer
	n = length(r)
	ln = Int(ceil(log(n)))
	supp = [rand(Uniform(-1,1))]
	exponents = [0:1:n-1;]
	ai = 2*sum((supp[1].^exponents).*r) - r[1]
	Bi = (1 + supp[1]^2)/(1-supp[1]^2)

	weight = [ai/Bi]

	conv = false
	count = 0
	while(count < 100 && !conv)
		validmeasure = false
		proposed = weight
		while(!validmeasure)
			B = zeros(length(supp), length(supp))
			for i in 1:length(supp)
				for j in 1:length(supp)
					B[i,j] = (1+supp[i]*supp[j])/(1-supp[i]*supp[j])	
				end
			end
			c = zeros(length(supp))
			for i in 1:length(supp)
				c[i] = 2*sum((supp[i].^exponents) .* r) - r[1]
			end

			prob = LinearProblem(B,c)
			sol = LinearSolve.solve(prob)
			unrestweight = sol.u
			
			if(all( >=(0), unrestweight))
				validmeasure = true
				weight = unrestweight
			else
				t = proposed./(-unrestweight .+ proposed)
				pop!(t)
				t[t .< 0] .= typemax(Int)
				t[t .> 0] .= typemax(Int)
				bd = findmin(t)
				proposed = (1-bd[1]).*proposed + bd[1] .* unrestweight
				deleteat!(supp, bd[2])
				deleteat!(proposed,bd[2])
			end

		end
		model = SOSModel(solver)
	@polyvar x
		
	f = 0
	for i in 1:length(supp)
		g=1
		for j in 1:length(supp)
			if( j != i)
				g = g* (1 - x*supp[j])
			end
		end
		f = f + weight[i]*(1 + x*supp[i])*g
	end
		
		
		h = f + (-2 * sum(r .* x.^exponents) + r[1])*prod((1 .- x.*supp))
		d = prod((1 .- x.*supp))
		
		print(h)
		print(d)
	 	S = @set x-a >= 0 && b-x >= 0
		@variable(model,s)
		
		@constraint(model,c, h >= s*d, domain = S)
		@objective(model, Max, s)
		optimize!(model)
		
		if(abs(objective_value(model)) < tol)
			conv = true
		end
		
		v = moment_matrix(model[:c])
		
		pt = atomic_measure(v,1e-4)
		print(pt)
		
		if(length(pt.atoms[1].center) == 1)
			append!(supp, pt.atoms[1].center)
			append!(pts, pt.atoms[1].center)
		else
			append!(supp, pt.atoms[1].center[2])
			append!(pts, pt.atoms[1].center[2])
		end
		append!(weight, 0)
		
		count = count + 1
	end
	
	return(supp, weight,pts)
end


SR = function(supp, weight,r)
		validmeasure = false
		current = weight
	exponent = [0:1:length(r)-1;]
		while(!validmeasure)
			B = zeros(length(supp), length(supp))
			for i in 1:length(supp)
				for j in 1:length(supp)
					B[i,j] = (1+supp[i]*supp[j])/(1-supp[i]*supp[j])	
				end
			end
			c = zeros(length(supp))
			for i in 1:length(supp)
				c[i] = 2*sum((supp[i].^exponent) .* r) - r[1]
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
				deleteat!(supp, bd[2])
				deleteat!(current,bd[2])
			end

		end
	return(supp, weight)
end


estimate_poly = function(i,r)
	m = Int(ceil(exp(1+1/exp(1))*log(10^12)))
	t = Int(floor(2^abs(i) * log(10^12)))
		a0 = (1- 2.0^-abs(i))
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



plot_poly = function(s)
	@polyvar x 
	exponents = [0:1:s[3]-1;]
	p = sum(s[1] .* (x - s[2]).^exponents)
	points = [0:0.01:0.5;]
	plot(points, p.(points))
end


estim_poly = function(s,x)
	exponents = [0:1:s[3]-1;]
	return sum(s[1] .* (x-s[2]).^exponents)
end




momentLSmod = function(r, delta,supp, weight, tol, graph = false)
	id = append!([1:1:13;], [-13:1:-1;])
	dictionary = Dict(id .=> [estimate_poly(i,r) for i = append!([1:1:13;],[-13:1:-1;])])
	pts = [-1+delta:0.01:1-delta;]
	n = length(r)
	exponents = [0:1:n-1;]
	conv = false
	count = 0
	while(count < 25 && !conv)
		SRstep = SR(supp, weight,r)
		supp = SRstep[1]
		weight = SRstep[2]
	
		
		points = grad_optimize(r, dictionary, supp, weight,delta)
		index = findmin(points[1])[2]
		if(findmin(points[1])[1] > -tol)
			conv = true
		end
		append!(supp, points[2][index])
		append!(weight, 0)
		if(graph == true)
            a = zeros(length(pts))
            b = zeros(length(pts))
            for i in 1:length(a)
                a[i] = -2*sum(r.* pts[i].^exponents) + r[1]
                b[i] = sum(weight.*(1 .+ pts[i].*supp)./(1 .- pts[i].*supp))
            end
            val = a+b
            if(count == 1)
                display(plot(pts, val))
            else
                display(plot!(pts,val))
            end
        end
            
		count = count + 1
	end
	
	return(supp, weight)
end



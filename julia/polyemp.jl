include("momentLS.jl")

solve_LDA = function(f,g)
    optim = zeros(26)
    id = append!([1:1:13;], [-13:1:-1;])
	p = Dict(id .=> [estimate_poly(i,f) for i = append!([1:1:13;],[-13:1:-1;])])
    time = 0
    i = 1
    for(ind in id)
        model = SOSModel(Clarabel.Optimizer)
        set_string_names_on_creation(model, false)
		@polyvar x
		fx = sum(p[ind][1] .* transform(x, ind) .^[0:1:(length(p[ind])-1);])
        gx = sum(g .* transform(x,ind) .^ [0,1:4;])

        S = @set x >= 0 && 1-x >= 0

        @variable(model,s)
			
		@constraint(model,c, fx >= s*gx, domain = S)
		@objective(model, Max, s)
		set_silent(model)
		optimize!(model)
        time += solve_time(model)
        optim[i] = objective_value(model)
    end
    min = findmin(optim)[1]
    return(time, min)
end

times = zeros(length(10:25:250;))
for i in [10:25:250;]

    @polyvar x
    f_c = rand(i)
    g_c = rand(5)
    f = sum(f_c .* x.^[1:1:i;])
    g = sum(g_c .* x.^[0:1:4;])
    S = @set x>=0 && x <= 1
    model = SOSModel(Clarabel.Optimizer)
    @variable(model, s)
    @objective(model, Max, s)
    @constraint(model, c, f -s*g>= 0,domain = S)
    optimize!(model)
    times[Int((i-10)/25) + 1] = solve_time(model)
end

println(times)
plot([10:25:250;], times)
savefig("timeplt.png")

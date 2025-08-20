using Random
import Clarabel
using SumOfSquares
using DynamicPolynomials
using Plots

times = zeros(length(10:25:250;))
for i in [10:25:250;]

    @polyvar x
    f = sum(rand(i) .* x.^[1:1:i;])
    g = sum(rand(5) .* x.^[0:1:4;])
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
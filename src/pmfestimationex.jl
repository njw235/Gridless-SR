include("SR1.jl")
using RCall
R"library(momentLS)"
R"measlist = vector('list', length = 50)"
Random.seed!(1234);
alphas = rand(10)
function sim_data(n, option)::AbstractArray{Integer}
    data = []
    for i in 1:n
        if(option == 1)
            r = rand(Uniform(0,1))
            if(r < 1/3)
                append!(data, rand(Geometric(0.8)))
            else
                append!(data, rand(Geometric(0.6)))
            end
        
        elseif(option == 2)
            r = rand(Uniform(0,1))
            if(r < 1/4)
                append!(data, rand(Geometric(0.9)))
            elseif(r > 1/4 && r < 3/4)
                append!(data, rand(Geometric(0.7)))
            else
                append!(data, rand(Geometric(0.2)))
            end
        
        elseif(option == 3)
            r = rand(Uniform(0,1))
            append!(data, rand(Geometric(1-r)))
        else
            r = rand(Uniform(0,1))
            for j in 1:10
                if r > (j-1)/10 && r < j/10
                    append!(data, rand(Geometric(alphas[j])))
                end
            end
        end
    end
    return(data)
end




# Simulation studies
# Using the two pmfs from the paper given
# Two finite mixtures of geometrics

i = parse(Int64, ARGS[1])

function pmft4(x)
    sum(0.1 .* alphas .* (1 .- alphas).^x)
end
function pmft3(x)
    1/((x+1)*(x+2))
end

function pmft2(x)
    0.25*0.9*0.1^x + 0.5*0.7*0.3^x + 0.25*0.2*0.8^x
end

function pmft1(x)
    (1/3)*0.8*0.2^x + (2/3)*0.6*0.4^x
end

N = 1000
errors = zeros(50)
oerrors = zeros(50)
for j in 1:50
    p = sim_data(N,i)
    r = Empirical(p)
    if(i == 1)
        d = 0.01
    elseif(i == 2)
        d = 0.01
    elseif(i==3)
        d = 0.001
    else
        d = 1- findmax(alphas)[1]
    end
    supp = [0.1]
    weight = [0.5]
    m = SR1_gridless(r, d,"weighted", "N", "LDA", 100, [0.0],[0.0], 1e-9)
    supp = m[1]
    weight = m[2]
    function pmf(x)
            sum(weight .* supp.^x)
    end        
     x = [0:1:10000000;]
    if(i == 1)
            errors[j] = sum((pmf.(x) .- pmft1.(x)).^2)               
    elseif(i == 2)
            errors[j] = sum((pmf.(x) .- pmft2.(x)).^2)               
    elseif(i==3)
            errors[j] = sum((pmf.(x) .- pmft3.(x)).^2)   
    else
            errors[j] = sum((pmf.(x) .- pmft4.(x)).^2) 
    end
    @rput j
    @rput supp
    @rput weight
    @rput r
    R"err = L2diff_L2Moment(r, supp, weight)"
    @rget err
    oerrors[j] = err
    println(errors[j], oerrors[j])
    R"measlist[[j]] = cbind(supp, weight)"
end

@rput i
R"saveRDS(measlist, paste(i, 'pmfmeas.rds', sep = ''))"


open(string("gridlesspmf",ARGS[1], ".txt"), "w") do io
    println(io, errors)
    println(io, oerrors)
end

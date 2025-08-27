include("mixingMeasure.jl")

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
        
        else
            r = rand(Uniform(0,1))
            append!(data, rand(Geometric(1-r)))
        end
    end
    return(data)
end




# Simulation studies
# Using the two pmfs from the paper given
# Two finite mixtures of geometrics





function pmft2(x)
    0.25*0.9*0.1^x + 0.5*0.7*0.3^x + 0.25*0.2*0.8^x
end

function pmft1(x)
    (1/3)*0.8*0.2^x + (2/3)*0.6*0.4^x
end
errordict = Dict()

Random.seed!(1234);

for N in [50,100,500,1000,5000,10000]
    errors = zeros(2,20)
    werrors = zeros(2,20)
    for j in 1:20
        for i in 1:2

            p = sim_data(N,i)
            r = Empirical(p)

            if(i == 1)
                d = 0.4
            elseif(i == 2)
                d = 0.1
            end

            supp = [0.1]
            weight = [0.5]

            m = mixingmeasure(r, d, supp, weight, 1e-8)

            supp = m[1]
            weight = m[2]
            function pmf(x)
                sum(weight .* supp.^x)
            end
            
            x = [0:1:10000;]

            if(i == 1)
                errors[i,j] = sum((pmf.(x) .- pmft1.(x)).^2)
               
            elseif(i == 2)
                errors[i,j] = sum((pmf.(x) .- pmft2.(x)).^2)
               
            else
                errors[i,j] = sum((pmf.(x) .- pmft3.(x)).^2)
                
            end
        end

    end

    errorlist = reduce(+, eachcol(errors)) ./ size(errors,2)

    
    errordict[N] = errorlist 

end



open("gridlesspmf.txt", "w") do io
    print(io, errordict)
end

using Random
using Distributions
using RCall

#Adding in the grid based implementation of the support reduction algorithm
# Documentation and code can be found at https://people.math.ethz.ch/~fadouab/ComparisonEstimSR.r
gs = ARGS[2]
@rput gs
R"library(momentLS)"
function sim_data(n, option)::AbstractArray{Integer}
  data = []
  for i in 1:n
    if (option == 1)
      r = rand(Uniform(0, 1))
      if (r < 1 / 3)
        append!(data, rand(Geometric(0.8)))
      else
        append!(data, rand(Geometric(0.6)))
      end

    elseif (option == 2)
      r = rand(Uniform(0, 1))
      if (r < 1 / 4)
        append!(data, rand(Geometric(0.9)))
      elseif (r > 1 / 4 && r < 3 / 4)
        append!(data, rand(Geometric(0.7)))
      else
        append!(data, rand(Geometric(0.2)))
      end

    else
      r = rand(Uniform(0, 1))
      append!(data, rand(Geometric(1 - r)))
    end
  end
  return (data)
end




# Simulation studies
# Using the three pmfs from the paper given
# Two finite mixtures of geometrics and a mixture model with
# Uniform distribution as the mixing measure 


function pmft3(x)
  1 / ((x + 1) * (x + 2))
end
function pmft2(x)
  0.25 * 0.9 * 0.1^x + 0.5 * 0.7 * 0.3^x + 0.25 * 0.2 * 0.8^x
end

function pmft1(x)
  (1 / 3) * 0.8 * 0.2^x + (2 / 3) * 0.6 * 0.4^x
end

i = ARGS[1]
Random.seed!(1234);
errors = zeros(50)
for j in 1:50
  p = sim_data(1000, i)

  if (i == 1)
    d = 0.4
  elseif (i == 2)
    d = 0.1
  else
    d = 0.002
  end

  supp = [0.1]
  weight = [0.5]
  @rput p
  R"counts=table(p)"
  R"pBar=rep(0,max(p)+1)"
  R"pBar[as.numeric(names(counts))+1]=counts/sum(counts)"
  R"alphaGrid=momentLS::makeGrid(upper_threshold = 1-1e-4,nX=as.integer(gs),cm=TRUE,scale='equidist')"
  R"XtX=1/(1-outer(alphaGrid,alphaGrid)) #use 1/(1-xy) inner product instead of (1+xy)/(1-xy)"
  R"s_alpha=sqrt(diag(XtX))"
  R"Xtr_0=Xtr_cpp(alphaGrid,pBar,standardization=FALSE,epsilon = 1e-10)"
  R"Xtr=pBar[1]/2+Xtr_0/2"
  R"XtX_standardized=diag(1/s_alpha)%*%XtX%*%diag(1/s_alpha)"
  R"Xtr_standardized=Xtr/s_alpha"
  R"fit=SR1(r=pBar,delta=NULL,alphaGrid = alphaGrid,
        precomputed=list(s_alpha=s_alpha,
                         XtX=XtX_standardized,
                         Xtr=Xtr_standardized,
                         alphaGrid=alphaGrid,
                         input=pBar))"
  R"supp = fit$support"
  R"weight = fit$weights/ (1-fit$support)"
  @rget supp
  @rget weight
  function pmf(x)
    sum((1 .- supp) .* weight .* supp .^ x)
  end

  x = [0:1:10000;]

  if (i == 1)
    errors[j] = sum((pmf.(x) .- pmft1.(x)) .^ 2)

  elseif (i == 2)
    errors[j] = sum((pmf.(x) .- pmft2.(x)) .^ 2)
  else
    errors[j] = sum((pmf.(x) .- pmft3.(x)) .^ 2)
  end
  println(errors[j])
end


open(string(gs,"gridpmf", i, ".txt"), "w") do io
  print(io, errors)
end

include("momentLS.jl")
R"library(momentLS)"
# Initialize a moment sequence supported by a single atom

Random.seed!(1234)
atom = rand(Uniform(0,.9))
weight = rand(Uniform(0,3))

n = 100

y = weight .* atom .^[0:1:n-1;]

dhat = 1-atom -0.05

m = momentLSmod(y, dhat, [0.0],[0.0], 1e-8)

@rput y

supp = m[1]
w= m[2]
pop!(supp)
pop!(w)
@rput supp
@rput w

R"diff = L2diff_L2Moment(y, supp, w)"

@rget diff

print("The estimated support is: ", supp, " The estimated weights are: ", w)
print("The true support is: ", atom, " The true weight is: ", weight)
print("The loss function value for our estimator is: ", diff)
using Random
using Distributions
using RCall

#Adding in the grid based implementation of the support reduction algorithm
# Documentation and code can be found at https://people.math.ethz.ch/~fadouab/ComparisonEstimSR.r

R"interval = seq(0,1,10^{-3})"
R"Emp <- function(x) {  
	  # compute the empirical frequence of x
	  n = length(x)
	  supp <- seq(0,max(x),by=1)
	  count <- rep(0,length(supp))
	  for(i in 1:length(count)) {
	    count[i] <- sum(x==supp[i])
	  }
	  
	  return(list(supp=supp,freq=count/n))
	  
	  
	}"
R"InitializeC = function(x, alpha){
  V  = 0:max(x)
  T = sum(Emp(x)$freq*alpha^V)
  
  return((1+alpha) * T)
  
}
"
R"
Interm = function(x,alpha){
  V=0:max(x)
  return(sum(Emp(x)$freq*alpha^V))
}"
R"DirDeriv = function(alpha, x, S, C){
  # C is the vector of weights  of the current iterate
  # S is the support vector of the current iterate
  
  T1 = sum(C* (1-S)/(1-S*alpha))
  T2 = Interm(x=x, alpha=alpha)
  
  return((1-alpha)* (T1 - T2))
  
}"
R"FindMinfun = function(x, S, C){
  
  V = 0:max(x)
  Out1 = rep(0, length(interval))
  Out2 = rep(0, length(interval))
  
  for(j in 1:length(S)){
    
    Out1 = Out1 + C[j]*(1-S[j])/(1-S[j]*interval)
  }
  
  
  for(j in 1:length(V)){
    
    Out2 = Out2 + Emp(x)$freq[j]*interval^{V[j]}
  }
  
  Out = (1-interval)*(Out1 - Out2)
  Ind =match(min(Out), Out)
  
  return(list(minimum=interval[Ind],objective=min(Out)))
  
}"
R"SolveUncons  = function(x, S){
  
  M = matrix(0, length(S), length(S))
  for(i in 1:length(S)){
    
    for(j in 1:length(S)){
      M[i,j]  = (1-S[j])/(1-S[i]*S[j])
    }
    
  }
  #M1=rep(1,length(S))%*%t(1-S)
  #M2 = 1-S%*%t(S)
  #M =M1/M2
  
  B= rep(0, length(S))
  for(i in 1:length(S)){
    B[i]  = Interm(x=x, alpha=S[i])
  }
  svdM = svd(M)
  uM = svdM$u
  vM = svdM$v
  diagM =diag(1/svdM$d)
  Sol1 = solve(uM, B)
  Sol2 = diagM%*%Sol1
  Sol  = vM%*%Sol2
  
  #Sol = solve(M, B)
  Sol= as.vector(Sol)
  List = list(Sol = Sol, Check1=as.vector(M%*%Sol)-B)
  
  return(List$Sol)
}"
	
R"ReduceSuppfun = function(S1, C1, S2, C2){
  
  S.m = c(S1, S2)
  S.m = unique(sort(S.m))
  
  # S1 will always be the support vector of the previous permissible iterate
  #  S2 will be the current support vector. An example:
  
  # S1=  c(0.2, 0.3, 0.4)
  # S2 = c(0.1, 0.2, 0.3, 0.4)  (this means that the point 0.1 was added to S1
  # Then S.m = c(0.1, 0.2, 0.3, 0.4) 
  
  S2.m = rep(10, length(S.m)) 
  
  # S2.m is a vector filled with the value 10. This choice is really arbitrary.
  
  C1.m = rep(0, length(S.m))  # we will store weights in this vector
  
  C2.m = rep(0, length(S.m))  # idem
  
  ind1 = pmatch(S1, S.m)
  ind2 = pmatch(S2, S.m)
  
  # The vector ind1 identifies the indices where S.m has the same values as in S1
  # The same thing for ind1
  # In the example given above we have: ind1 = c(2,3,4) and ind2 = c(1,2,3,4)
  # 
  
  
  
  C1.m[ind1]  = C1   # we fill the right positions with the weights from C1
  C2.m[ind2]  = C2   # idem
  
  S2.m[ind2]  = S2    # we fill the right positions with the support points from S2
  
  Lambda = C1.m/(C1.m-C2.m)  # we compute this ratio (componentwise) as the mimimal value of the positive
  # elements of Lambda gives the right convex combination (which will delete
  # one of the support point of S2
  
  index = match(min(Lambda[Lambda > 0]), Lambda)  # We look for the position of the minimum
  
  S2new  = S2.m[-index]   # We remove that support points from the augmented support vector corresponding to S2
  S2new = S2new[S2new < 10]  # We get rid of the values 10 to get finally the reduced suppport vector 
  
  return(S2new)
}"
R"CompMonLSE= function(x,alpha0,epsilon){
  
  S0=alpha0   # any intial value of choice
  
  C0 = InitializeC(x=x, alpha=alpha0)  # the initial optimal weight
  #print(C0)
  
  Resmin = FindMinfun(x=x, S=S0,C=C0)  # Looking for the optimal Dirac direction in which we will go 
  thetamin=Resmin$minimum              # This is the new support point that will be added
  valmin = Resmin$objective
  
  count =0
  Count = 0
  while(abs(valmin) > epsilon ){   # a stopping condition
    count <- count + 1	 
    Snew=c(S0,thetamin)  # The new suport
    Snew = sort(Snew)
    
    
    Cnew=SolveUncons(x=x,S=Snew)  # Solve the unconstrained LS problem
    
    min.C <- min(Cnew) 
    #Count <- 0
    while(min.C < 0){  # condition indicating that the unconstrained solution is not permissible
      # and that we should go into the reduction step
      
      Count <- Count+1

      Snew = ReduceSuppfun(S1=S0, C1=C0, S2=Snew, C2=Cnew)
      
      if(length(Snew) > 1){  # two expressions for the solutions depending on whether we are left with 
        # more than 1 point or just one
        Cnew = SolveUncons(x=x, S=Snew)
      }
      else{
        Cnew=InitializeC(x=x, alpha=as.numeric(Snew))
        
      }  #  the program will go out of the inner loop (reduction loop) if 
      # and only if the solution is permissible
      
      min.C = min(Cnew)
      
      
    }#min.C
    
    S0= Snew  # uppdate the support vector
    C0= Cnew  # updfate the weight vector
    
    Resmin = FindMinfun(x=x, S=S0,C=C0)  # compute the new minimizer of the directional derivative
    # to build up a new support or to stop if the stopping 
    #crierion is met. 
    
    thetamin=Resmin$minimum
    valmin = Resmin$objective
    #print(valmin)
  }#valmin
  
  return(list(S0=S0,C0=C0, valmin=valmin))
  
}"



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
	# Using the three pmfs from the paper given
	# Two finite mixtures of geometrics and a mixture model with
	# Uniform distribution as the mixing measure 
	
	
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
	
	            if(i == 1)
	                d = 0.4
	            elseif(i == 2)
	                d = 0.1
	            else
	                d = 0.0002
	            end
	
	            supp = [0.1]
	            weight = [0.5]
				@rput p
	            R"m = CompMonLSE(p, 0.1, 1e-8)"
	
	            R"supp = m$S0"
	            R"weight = m$C0"
				@rget supp
				@rget weight
	            function pmf(x)
	                sum((1 .- supp) .* weight .* supp.^x)
	            end
	            
	            x = [0:1:10000;]
	
	            if(i == 1)
	                errors[i,j] = sum((pmf.(x) .- pmft1.(x)).^2)
	               
	            elseif(i == 2)
	                errors[i,j] = sum((pmf.(x) .- pmft2.(x)).^2)
	        end
	
	    end
	
	    errorlist = reduce(+, eachcol(errors)) ./ size(errors,2)
	
	    
	    errordict[N] = errorlist 
	
	end
	
	print(errordict)
end

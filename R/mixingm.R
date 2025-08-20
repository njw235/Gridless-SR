interval = seq(0,1,10^{-5})

#--------------------------------------------------------------
#This function empirical distribution based on 
# observed data "x". It returns a vector of length equal 
# to max(x) +1 because we deal here with data of the form 
# 0,1,..., max(x). If some value is not observed, the function
#returns 0.
#--------------------------------------------------------------

Emp <- function(x) {  
  # compute the empirical frequence of x
  n = length(x)
  supp <- seq(0,max(x),by=1)
  count <- rep(0,length(supp))
  for(i in 1:length(count)) {
    count[i] <- sum(x==supp[i])
  }
  
  return(list(supp=supp,freq=count/n))
  
  
}

#----------------------------------------------------------
#For given observed data "x" and given candidate "alpha" in
#[0,1] (or rather the grid above),  the following function
# gives the corresponding optimal weight: the function
#returns "\hat{pi}" such that the completely monotone
# function "\hat{\pi} (1-\alpha) \alpha^i, i \in \NN", is the
#the LS projection of the empirical estimator obtained by
# "Emp" 
#---------------------------------------------------------

InitializeC = function(x, alpha){
  V  = 0:max(x)
  T = sum(Emp(x)$freq*alpha^V)
  
  return((1+alpha) * T)
  
}


##########This intermediate/auxiliary function is used for computing the directional
##########derivative below


Interm = function(x,alpha){
  V=0:max(x)
  return(sum(Emp(x)$freq*alpha^V))
}


######Now, given a current candidate for the Least Squares estimtion problem
#########with support vector "S" and weights vector "C", the following function computes the directional
#########derivative of the LS criterion in the direction of a geoemtric pmf with prob "\alpha"
############Note that this function calls this "Interm" function to make things less cumbersome.

DirDeriv = function(alpha, x, S, C){
  # C is the vector of weights  of the current iterate
  # S is the support vector of the current iterate
  
  T1 = sum(C* (1-S)/(1-S*alpha))
  T2 = Interm(x=x, alpha=alpha)
  
  return((1-alpha)* (T1 - T2))
  
}

}


##########We need to find the minimum in "\alpha"  of the directional derivative
##########defined above. That's what this function is doing. I have decide to look
##########for the minimum using my hands and not giving this task to some optimizing  
#############function

FindMinfun = function(x, S, C){
  
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
  
}


##########The following function solves the LS problem in the "unconstrained" set, that is, 
##########for a given support vector "S", it finds the optimal weight vector "C" so that the
########## elements of "C" are allowed to be negative. The system to solve is linear: M C = B
##########where "M" is a matric which can be very ill-conditionned. Here, WE NEED HIGH PRECISION.
###########The method I use to solve the system can be changed of course. The only things that should
###########REMAIN are "M" and "B".


SolveUncons  = function(x, S){
  
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
}



################The following function makes the support reduction step works. When the unconstrained
##############solution is not permissible, that is, some elements of "C" are negative, the function
############## tells us which support point to remove from "S" by computing a special convex combination
############## of the previous permissible iterate and the new non-permissible one. 
##################The way this function works is described in comments I've put in the function itself

ReduceSuppfun = function(S1, C1, S2, C2){
  
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
}


####This function is the main one which should give the completely monotone LSE


#-------------------------------------------------------------------------
#The Support reduction algorithm for computing the completely monotone LSE
# if epsilon > 10^{-9}, then one does not get stuck. So, typically one can
# choose epsilon no smaller than 10^{-9} but no bigger than 10^{-7}
#alpha0 can be given also to the user, the default is 0.1 
#-------------------------------------------------------------------------

CompMonLSE= function(x,alpha0,epsilon){
  
  S0=alpha0   # any intial value of choice
  
  C0 = InitializeC(x=x, alpha=alpha0)  # the initial optimal weight
  #print(C0)
  
  Resmin = FindMinfun(x=x, S=S0,C=C0)  # Looking for the optimal Dirac direction in which we will go 
  thetamin=Resmin$minimum              # This is the new support point that will be added to "\alpha0" 
  valmin = Resmin$objective
  
  count =0
  Count = 0
  while(abs(valmin) > epsilon ){   # a stopping condition
    count <- count + 1	
    #cat("Main loup numb = ",count,"\n")  
    Snew=c(S0,thetamin)  # The new suport
    Snew = sort(Snew)
    
    
    Cnew=SolveUncons(x=x,S=Snew)  # Solve the unconstrained LS problem
    
    min.C <- min(Cnew) 
    #Count <- 0
    while(min.C < 0){  # condition indicating that the unconstrained solution is not permissible
      # and that we should go into the reduction step
      
      Count <- Count+1
      #cat("Sub loup numb = ",Count," of the main loop numb=", count, "\n")
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
  
}
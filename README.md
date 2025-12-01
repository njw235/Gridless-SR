# A Gridless Support Reduction Algorithm For Finite Mixture Problems

This repository contains the code and experiments for the paper "A Gridless Support Reduction Algorithm For Finite Mixture Problems".

# Setting Up Environment

1. Clone the repository:
```
$ git clone https://github.com/njw235/Gridless-SR.git
$ cd Gridless-SR
$ julia
```
2. Activate the project environment for the required dependencies:
```
> import Pkg
> Pkg.activate(".")
> Pkg.instantiate()
```

# File Content

This project consists of R files and julia files all in the src folder. The relevant files for replicating the tables in the paper are: 
1. Blassosim.jl - contains code to run the simulations of the Blasso example
2. MHsim.jl - contains code to run simulations of Metropolis-Hastings examples
3. Testsmoments.jl - Contains code to run the AR(1) simulations.
4. pmfestimationex.jl - contains code to run the pmf estimation simulations.
5. time_comparison.R - contains code comparing the runtime LDA to no approximation for the grid based SR. 

# Running a given file

## Running a julia file

Once you have activated the project environment and have the required dependencies from the file, you can run any of the given julia files from the command line:

```
$julia julia/filename.jl
```

Note that in the case of the simulation examples, all except the Metropolis-Hastings example expect an input in the command line to know which chain/sample the estimation procedure should be run on.
Specifically Blassosim.jl expects a number between 1-9, Testmoments.jl 1-5, pmfestimationex.jl 1-4. 
## Running a pluto notebook

To run the pluto notebooks execute the following code:

```
>using Pkg
>Pkg.add("Pluto")
>using Pluto
>Pluto.run()
```

This will open a notebook environment. From there, go to the "Open a notebook" path and input the filepath in which you have the local version of the pluto notebook you wish to run.
From there, just run each cell individually and see the output. 


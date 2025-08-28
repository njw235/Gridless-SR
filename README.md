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

This project consists of R files, raw julia files, and Pluto Notebooks.

R files:

1. time_comparison.R - This file contains the implementation for comparing the runtime and performance of the grid version of the support reduction algorithm for AR(1) chains with and without the low-degree polynomial approximation

Julia files:

1. momentLS.jl: This file contains the implementation of the low-degree polynomial approximation and the gridless support reduction algorithm.
2. mixingMeasure.jl: This file contains the implementation as above but for the pmf estimation problem.
3. polyemp.jl: This file contains the empirical testing of the time complexity of the SOS problem described in section 5.1.
4. pmfestimationex: This file contains the empirical testing for pmf estimation in section 5.2.
5. auxMCMCtest.jl: This file contains the empirical testing for the BLASSO and MH chains in section 5.3.
6. Testmomentls.jl: This file contains the empirical testing for the AR(1) example in section 5.3.
7. SmallExample.jl: This file contains a small example showing the fitting of a moment sequence that is easily runnable as opposed to the longer run times of the other files.

Pluto Notebooks:

1. AutocovariancePlot.jl: This file contains a runnable pluto notebook to reproduce the autocovariance plots.
2. pmfplots.jl: This file contains a runnable pluto notebook to reproduce the pmf error plots.
3. RunnableExample.jl: This file contains a small runnable example from SmallExample.jl but in the form of a pluto notebook for ease of readibility.

# Running a given file

Once you have activated the project environment and have the required dependencies from the file, you can run any of the given julia files from the command line:

```
$julia julia/filename.jl
```

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


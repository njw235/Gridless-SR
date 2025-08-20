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

This project consists of both R and julia files

Julia files:

1. momentLS.jl: This file contains the implementation of the low-degree polynomial approximation and the gridless support reduction algorithm.
2. mixingMeasure.jl: This file contains the implementation as above but for the pmf estimation problem.
3. polyemp.jl: This file contains the empirical testing of the time complexity of the SOS problem described in section 5.1.
4. pmfestimationex: This file contains the empirical testing for pmf estimation in section 5.2.
5. auxMCMCtest.jl: This file contains the empirical testing for the BLASSO and MH chains in section 5.3.
6. Testmomentls.jl: This file contains the empirical testing for the AR(1) example in section 5.3.
7. SmallExample.jl: This file contains a small example showing the fitting of a moment sequence that is easily runnable as opposed to the longer run times of the other files.

# Running a given file

Once you have activated the project environment and have the required dependencies from the file, you can run any of the given files from the command line:

```
$julia julia/filename.jl
```


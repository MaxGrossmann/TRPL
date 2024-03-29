# Time-Resolved Photoluminescence (TRPL)
**Generalized modelling of photoluminescence transients**

If you use our data or software as part of your research, teaching, or other activities, \
we would be grateful if you could cite our work.

Generalized Modeling of Photoluminescence Transients - https://doi.org/10.1002/pssb.202200339

## Disclaimer
The provided optimization scripts 
*trpl_opt_real.jl* and *trpl_opt_syn.jl* should only be used testing purposes
(output goes to "\TRPL\Optimizations\Test"). For a real optimziation one would start multiple optimization runs 
using an additional script by changing lines 12-19 to ARGS[1], ARGS[2],... in *trpl_opt_real.jl* and *trpl_opt_syn.jl*

When first using the code you need to instantiate the packages using Pkg.instantiate() inside the Code directory.
(for more information visit https://docs.julialang.org/en/v1/). 

For best performance
one would do one optimization run per script distributed over multiple CPU cores (using 1-2 cores per optimization run), saving
the results of each run in a seperate directory. We obtained our results in this way by submitting multiple batch jobs on our local computer cluster. An example can be found in the following directory: "\TRPL\Optimizations\Real\sample1".

## Code
contains:

- main optimization code (trpl_opt_real.jl and trpl_opt_syn.jl)
- analysis of optimization results (trpl_opt_analysis_real.jl & trpl_opt_analysis_syn.jl)
- index of reflexion interpolation (ior_interp.jl)
- calculation of transmission coefficents using transfer-matricies (transm_calc.jl)
- calculation of experimental initial conditions (eta0_calc.jl)
- mono-exponential fits (mono_exp_fits.jl)
- generation of synthetic data sets (syn_data_gen.jl)
- Manifest.toml and Project.toml include information about dependencies, package names, etc.

## Data

contains:

- subfolders Real & Syn
    - Real: data from sample 1 and sample 2
    - Syn: generate synthetic data sets (subfolder CSV) and parameters used for there creation (subfolder Params)

structure of real data sets:

- colums: time, transient 1, transient 2, ..., transient N

- lines 1-3 (skipping the 1. column in each line) contain the following experimental parameters: 
    - laser power P
    - ND filter values
    - average background noise b

- line 4-end contains the binned experimental data (time, transient 1, transient 2, ..., transient N)

structure of synthetic data sets:

- colums: time, transient 1, transient 2, ..., transient N

- lines 1-3 (skipping the 1. column in each line) contain the following parameters: 
    - C scaling factor
    - ND filter values
    - background noise b

- line 4-end contains the binned synthetic data (time, transient 1, transient 2, ..., transient N)

## Optimizations

contains: 

- subfolders Real & Syn
    - Real: optimization results for sample 1 and 2 
    - Syn: optimization results for the synthetic data sets 1-5

Each folder contains 25 subfolders for each optimization run. In each of these folders
contains a .csv file with the optimization fitness, optimizied parameters, parameter bounds and the random seed.

## Results

Contains final results, parameter sets for all optimization runs, the data used and the fit curves as one .xlsx files for each real or synthetic data set.
The residual analysis can be found in separated .csv files. For the real data sets we included the results of the approximation of the experimental initial conditions as well as the results of the mono-exponential fits of the lowest transients of each real data set. 

In the results tab of the .xlsx-files one finds p_fit and p_opt. p_fit is the best parameter set found in all global optimizations.
p_opt is the final result obtained through a reparameterized local optimization starting at p_fit.

## Licencing 

Copyright (c) 2022 Max Großmann

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

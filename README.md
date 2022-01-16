# TRPL
**Beyond curve fitting: Modelling photoluminescence kinetics**

Code, Data and Results from the following publication: "insert DOI"

## Disclaimer
The data from our two samples was provided by the *Fraunhofer ISE*. The measurements
were made by *Klaus Schwarzburg* at the *Helmholtz-Zentrum Berlin*. The provided optimization scripts 
*trpl_opt_real.jl* and *trpl_opt_syn.jl* can be used in this repository for testing purposes only 
(output goes to "\TRPL\Optimizations\Test"). For a real optimziation one would start multiple optimization runs 
using an additional script by changing lines 12-19 to ARGS[1], ARGS[2],... in *trpl_opt_real.jl* and *trpl_opt_syn.jl*

When first using the code you need to instantiate the packages using Pkg.instantiate() inside the Code directory.
(for more information visit https://docs.julialang.org/en/v1/). 

For best performance
one would do one optimization run per script distributed over multiple CPU cores (using 1-2 cores per optimization run), saving
the results of each run in a seperate directory. We obtained our results in the described way by submitting multiple batch jobs on our local computer cluster. An example can be found in the following directory: "\TRPL\Optimizations\Real\sample1".

## Code
contains:

- main optimization code (trpl_opt.jl)
- analysis of optimization results (trpl_opt_analysis_real.jl & trpl_opt_analysis_syn.jl)
- index of reflexion interpolation (ior_interp.jl)
- calculation of transmission coefficents using transfer-matricies (transm_calc.jl)
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

- lines 1-3 (skipping the 1. colum in each line) contain the following experimental parameters: 
    - laser power P
    - ND filter values
    - average background noise b

- line 4-end contains the binned experimental data (time, transient 1, transient 2, ..., transient N)

structure of synthetic data sets:

- colums: time, transient 1, transient 2, ..., transient N

- lines 1-3 (skipping the 1. colum in each line) contain the following parameters: 
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

Contains final results, all optimization runs and the data used as .xlsx files for real and synthetic optimizations. 

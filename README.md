# TRPL
**Beyond curve fitting: Modelling photoluminescence kinetics**

Code, Data and Results from the following publication: "insert DOI"

## Disclaimer
The data from our two samples was provided by the *Fraunhofer ISE*. The measurements
were made by *Klaus Schwarzburg* at the *Helmholtz-Zentrum Berlin*. The files for sample 1 are called 
C2704_hzb..., for sample 2 C2702_hzb. The provided optimization scripts *trpl_opt_real.jl* and *trpl_opt_syn.jl*
can be used in this repository for testing purposes only (output goes to "\TRPL\Optimizations\Test"). For a real optimziation
one would start multiple optimization runs using an additional script by changing lines 12-19 to ARGS[1], ARGS[2],... 
(see https://docs.julialang.org/en/v1/manual/getting-started/) in *trpl_opt_real.jl* and *trpl_opt_syn.jl*. For best performance
one would do one optimization run per script distributed over multiple CPU cores (using 1-2 cores per optimization run), saving
the results of each run in a seperate directory. Using this procedure we obtained our results. An example can be found in the following directory: "\TRPL\Optimizations\Real\C2702_hzb".


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
    - Real: data from sample 1 (C2704) and sample 2 (C2702)
    - Syn: generate synthetic data sets (subfolder CSV) and parameters used for there creation (subfolder Params)

## Optimizations

contains: 

- subfolders Real & Syn
    - Real: optimization results for sample 1 and 2 
    - Syn: optimization results for the synthetic data sets 1-5

Each folder contains 25 subfolders for each optimization run. In each of these folders
are .log and .err files, a .lsf script to submit a batchjob at our local computing cluster and 
a .csv containing the optimization fitness, optimizied parameters, parameter bounds and the random seed.

## Results

Contains final results, all optimization runs and the data used as.xlsx files for real and synthetic optimizations. 

## Annotations

- PRNG = Mersenne Twister: for optimizations and Poisson Noise for synthetic data sets
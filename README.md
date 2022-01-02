# TRPL
Generalized Analysis of Time-resolved Photoluminescence Measurements

Code, Data and Results from the following publication: "insert DOI"

## Code
contains:

- main optimization code (trpl_opt.jl)
- analysis of optimization results (trpl_opt_analysis_real.jl & trpl_opt_analysis_syn.jl)
- index of reflexion interpolation (ior_interp.jl)
- calculation of transmission coefficents using transfer-matricies (transm_calc.jl)
- mono-exponential fits (mono_exp_fits.jl)
- generation of synthetic data sets (syn_data_gen.jl)

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

(PRNG = Mersenne Twister: for optimizations and Poisson Noise for synthetic data sets)
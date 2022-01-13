## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")

## physical constants
h = 6.62607015e-34
c = 299792458

## experimental parameters
rep_rate = 200e3
P        = [6.89e-8, 6.21e-4] # min P and max P of our laser
lambda   = 520*1e-9
fwhm     = 205*1e-6
d        = 2000*1e-9
T        = 0.4292 # calculated in transm_calc.jl
N_DA     = 5e16 # sample 1: N_DA = 5e16, sample 2: N_DA = 1e17

## main calculations
E        = P./rep_rate
N        = T*E/(h*c)*lambda
V        = pi*fwhm^2*d/4
n        = N/V * 1e-6
fwhm_int = 0.754 # from fwhm_int.jl script
eta_0    = round.(fwhm_int*n/N_DA,sigdigits=3)

## print result
println("\nη₀ ∈ [$(eta_0[1]),$(eta_0[end])]")

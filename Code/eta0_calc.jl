## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")
## constants
h = 6.62607015*1e-34
c = 299792458
## parameters
rep_rate = 200e3
P = [6.89e-8, 6.21e-4]
lambda = 520*1e-9
fwhm = 205*1e-6
d = 2000*1e-9
T = 0.4292
N_DA = 5e16
## calcuations
E = P./rep_rate
N = T*E/(h*c)*lambda
V = pi*fwhm^2*d/4
n = N/V * 1e-6
fwhm_int = 0.76 # from fwhm_int.jl script
eta_0 = round.(fwhm_int*n/N_DA,sigdigits=3)
## print result
println("\nη₀ ∈ [$(eta_0[1]),$(eta_0[end])]")

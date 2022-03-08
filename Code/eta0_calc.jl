## path to TRPL folder --> set it before running the code!
path_to_trpl = "D:\\Projects\\TRPL" # some path to the TRPL folder (just an example here)

## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")

## loading packages
using CSV, DataFrames

## directory
cd(joinpath(path_to_trpl,"Data\\Real"))

## output results? (true or false)
out = false

## parameters measurement data
name_measurement = "sample1"

if name_measurement == "sample1"
    N_DA = 5e16 # sample1: N_DA = 5e16
elseif name_measurement == "sample2"
    N_DA = 1e17 # sample2: N_DA = 1e17
end

## load measurement data
data = Matrix(CSV.read(name_measurement*".csv",DataFrame,skipto=2))

## laser power P 
P = round.(data[1,2:end],sigdigits=3)

## transmission coefficient (calculated in transm_calc.jl)
T = 0.4292 

## other experimental parameters
rep_rate = 200e3
lambda   = 520*1e-9
fwhm     = 205*1e-6
d        = 2000*1e-9

## physical constants
h = 6.62607015e-34
c = 299792458

## main calculations
E        = P./rep_rate
N        = T*E/(h*c)*lambda
V        = pi*fwhm^2*d/4
n        = N/V * 1e-6
fwhm_int = 0.754 # from fwhm_int.jl script
eta_0    = round.(fwhm_int*n/N_DA,sigdigits=3)

## print result
println("P = $(P)")
println("η₀ = $(eta_0)")


## write output
if out
    cd(joinpath(path_to_trpl,"Results\\Real"))
    file_name = "eta0_$(name_measurement).txt"
    if isfile(file_name)
        rm(file_name)
    end
    open(file_name,"w") do io
        println(io,"\nP = $(P)")
        println(io,"\nη0 = $(eta_0)")
    end
end
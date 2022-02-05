## path to TRPL folder
path_to_trpl = "D:\\Projects\\TRPL"

## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")

## loading packages
using XLSX, Distributions, Measurements

## directory
cd(joinpath(path_to_trpl,"Results\\Real"))

## name of real data set (sample 1: C2704, sample 2:C2702)
name = "C2702_hzb.xlsx"

## get xlsx and sort them
files =  filter(x->endswith(x, ".xlsx"), readdir())

## get results sheet
println(files[files .== name][1])
println("")
xf = XLSX.readxlsx(files[files .== name][1])
sh = xf["Results"]
num_pars = 10

## confidence level (1-α)
alpha = 0.05 # corresponds to a 95% confidence level
q = quantile(Normal(),1-alpha/2)

## Blakemore GaAs (valence band) 
N_V = 9.51e18

## get parameters
data = Array{Any}(undef,2,num_pars + 1)
data[1:2,2:end] = sh["C3:L4"]
data[2:end] ./= q
p = data[1,2:end] .± data[2,2:end]

## calculate band offset
offset = -25.7e-3 * log((p[4] * p[7]) / (p[5] * N_V))

## print output
println("band offset of shallow trapping state = $offset")
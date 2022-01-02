## path to TRPL folder --> set it before running the code!
path_to_trpl = "C:\\TRPL" # some path to the TRPL folder
## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")
## loading packages
using DataFrames, CSV, XLSX, LinearAlgebra, ForwardDiff, NLopt, Measurements, PGFPlots, ColorSchemes
## directory
cd(string(path_to_trpl,"\\Data\\Real\\CSV"))
## output results? (true or false)
out = false
## name of measurement 
name_measurement = "C2704_hzb"
## load measurement data
data_measurement = Matrix(CSV.read(string(name_measurement,".csv"),DataFrame,skipto=2))
## function for data processing
function data_processing(data)
    t = data[4:end,1]
    idx_end = sum(t .< 200)
    num_traj = size(data,2)-1
    dims = 9 + num_traj
    P = round.(data[1,2:end],sigdigits=3)
    ND = 10 .^ data[2,2:end]
    background = data[3,2:end]
    IPL = data[4:end,2:num_traj+1]
    tspan = (t[1], t[idx_end])
    N = length(t)
    idx_IPL = zeros(Bool,size(IPL))
    for i = 1:num_traj
	    idx_IPL[:,i] = IPL[:,i] .> 0
    end
    C_scaling = maximum(IPL)
    return t[1:idx_end], tspan, num_traj, dims, IPL[1:idx_end,:], idx_IPL[1:idx_end,:], P, ND, background, N, C_scaling 
end
## loading of experimental data
t, tspan, num_traj, dims, IPL, idx_IPL, P, ND, background, N, C_scaling = data_processing(data_measurement)
IPL = IPL[:,1]
background = background[1]
## mono-exponential model
model(p) = @. p[1] * exp(-t/p[2]) + background
## mle loss function 
function loss(par)
    sim = model(par)
    idx_used = (sim .> 0) .& (IPL .> 0)
    l = sum(sim[idx_used] .- IPL[idx_used] .* log.(sim[idx_used]))
    return l
end
## definition nlopt loss function
function loss_nlopt(par::Vector,grad::Vector)
    if length(grad) > 0
	    grad[:] = ForwardDiff.gradient(loss,par)
    end
    println("l = $(loss(par))")
    return loss(par)
end
## LBFGS (local optmization)
opt = Opt(:LD_LBFGS,2)
opt.min_objective = loss_nlopt
opt.ftol_abs = 1e-14
opt.maxtime = 600
(optf,p_opt,ret) = optimize(opt,[1000, 10])
## calculate standard errors
H_opt = ForwardDiff.hessian(loss,p_opt)
SE_opt = sqrt.(diag(inv(H_opt)))
## print results
println("\nτ_mono = $(round(p_opt[2],sigdigits=3)) ± $(round(SE_opt[2],sigdigits=3)) ns")
## get fit parameters with standard erros
cd(string(path_to_trpl,"\\Results\\Real"))
files = filter(x->endswith(x, ".xlsx"), readdir())
for i = 1:length(files)
    if occursin(Regex("$name_measurement"),files[i])
        global idx = i
    end
end
xf = XLSX.readxlsx(files[idx])
sh = xf["Results"]
p_fit = Float64.(sh["C3:K3"])
SE_fit = Float64.(sh["C4:K4"])
p_fit = p_fit .± SE_fit
## compair to effective lifetime 
tau_eff = round(sum([1/p_fit[i] for i = 1:6]) .^ -1,sigdigits=3)
println("\nτ_eff  = $tau_eff ns")
## write output
if out 
    file_name = "mono_exp_fit_$(name_measurement).txt"
    if isfile(file_name)
        rm(file_name)
    end
    open(file_name,"w") do io
        println(io,"\nτ_mono = $(round(p_opt[2],sigdigits=3)) ± $(round(SE_opt[2],sigdigits=3)) ns")
        println(io,"\nτ_eff  = $tau_eff ns")
    end
end
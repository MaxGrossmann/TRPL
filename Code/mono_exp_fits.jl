## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")
## loading packages
using DataFrames, CSV, XLSX, LinearAlgebra, ForwardDiff, NLopt, Measurements, PGFPlots, ColorSchemes
## directory
cd(string(homedir(),"\\OneDrive\\Uni_Master\\TRPL\\Data\\HZB\\CSV"))
## parameters measurement data
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
## mle cost function 
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
## LBFGS
opt = Opt(:LD_LBFGS,2)
opt.min_objective = loss_nlopt
opt.ftol_abs = 1e-14
opt.maxtime = 600
(optf,p_opt,ret) = optimize(opt,[1000, 10])
## calculate standard errors
H_opt = ForwardDiff.hessian(loss,p_opt)
SE_opt = sqrt.(diag(inv(H_opt)))
## plot data and fit 
IPL_fit = model(p_opt)
ylim = (10^-0.5,10^(ceil(log10(maximum(IPL)))))
plot_publication = PGFPlots.Axis(style="width=8.5cm, height=8.5cm, grid=major", ymode="log",
								 xlabel="t (ns)", ylabel="Counts", xmin=0, xmax=t[end],
								 ymin=ylim[1], ymax=ylim[2])
push!(plot_publication, Plots.Linear(t,IPL.+1e-3, style="red, smooth, solid", 
                                     mark="none", legendentry=L"$P_1$"))
push!(plot_publication, Plots.Linear(t,IPL_fit.+1e-3,style="black, smooth, solid",mark="none"))
file_name_tex= string(homedir(),"\\OneDrive\\Uni_Master\\TRPL\\Data\\Final\\Supplement\\Tikz\\mono_exp_fit_$(name_measurement).tex")
save(file_name_tex,plot_publication)
## print results
println("\nτ_mono = $(round(p_opt[2],sigdigits=3)) ± $(round(SE_opt[2],sigdigits=3)) ns")
## get fit parameters with standard erros
cd(string(homedir(),"\\OneDrive\\Uni_Master\\TRPL\\Data\\Final\\Real\\data"))
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
file_name = "mono_exp_fit_$(name_measurement).txt"
if isfile(file_name)
    rm(file_name)
end
open(file_name,"w") do io
    println(io,"\nτ_mono = $(round(p_opt[2],sigdigits=3)) ± $(round(SE_opt[2],sigdigits=3)) ns")
    println(io,"\nτ_eff  = $tau_eff ns")
end
## path to TRPL folder --> set it before running the code!
path_to_trpl = "D:\\Projects\\TRPL" # some path to the TRPL folder (just an example here)

## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")

## loading Packages
using OrdinaryDiffEq, DiffEqSensitivity, ForwardDiff, LinearAlgebra, 
      Random, Printf, Distributions, BlackBoxOptim, NLopt, DataFrames, CSV

## input parameters
# Parameters here are setup for a quick test! 
# For a real optimization set n_iter >= 2500000, max_time_loc_opt = 3600 to be save!
job_dir = joinpath(path_to_trpl,"Optimizations\\Test\\Job1") # job directory (just for tests)
name_measurement = "sample1" # name of measurement
n_iter = 2500 # number of optimization iteration and runs
n_runs =  1 # number of optimization runs
ftol_loc_opt = 1e-14 # tolerance for local optimization
max_time_loc_opt = 15 # maximum time in s for local optimization
#          τ_r0,  τ_nr0, τ_2t, τ_1t, τ_1e, τ_1d,    α,    β,    r,   C 
lb = vcat([   1,      1,    1,    1,    1,    1,    1,    1,    1,   1], 
           repeat([1e-6],10)) # lower parameter bounds
#          τ_r0,  τ_nr0, τ_2t, τ_1t, τ_1e, τ_1d,    α,    β,    r,   C 
ub = vcat([5000,   5000, 1000, 1000,  1e9, 1000,  1e5,  1e5, 1000, 1e6], 
           repeat([ 1.0],10)) # upper parameter bounds

## data directory
cd(joinpath(path_to_trpl,"Data\\Real"))

## load measurement data
data_measurement = Matrix(CSV.read(name_measurement*".csv",DataFrame,skipto=2))

## job directory
cd(job_dir)

## function for data processing (data structre is explained in README.md)
function data_processing(data)
    t = data[4:end,1]
    num_traj = size(data,2)-1
    num_pars = 10
    dims = num_pars + num_traj
    P = round.(data[1,2:end],sigdigits=3)
    ND = 10 .^ data[2,2:end]
    background = data[3,2:end]
    IPL = data[4:end,2:num_traj+1]
    tspan = (t[1], t[end])
    N = length(t)
    idx_IPL = zeros(Bool,size(IPL))
    for i = 1:num_traj
	idx_IPL[:,i] = IPL[:,i] .> 0
    end
    C_scaling = maximum(IPL)
    return t, tspan, num_traj, num_pars, dims, IPL, idx_IPL, P, ND, background, N, C_scaling 
end

## process measurement data
t, tspan, num_traj, num_pars, dims, IPL, idx_IPL, P, ND, background, N, C_scaling = data_processing(data_measurement)

## lower and upper bound dimension
lb = lb[1:dims]
ub = ub[1:dims]

## definition rate equations
function rate_equations(du, u, p, t)
    du[1] = -((1 - u[2]) / p[3]) * u[1] - ((1 - u[3]) / p[4]) * u[1] + (1 / p[5]) * u[3] * 1 / p[7] -
            ((1 + u[1]) / p[1]) * u[1] - ((1 + u[1]) / (p[2] * (1 + p[9] * u[1]))) * u[1]
    du[2] = ((1 - u[2]) / (p[3])) * p[8] * u[1]
    du[3] = ((1 - u[3]) / p[4]) * p[7] * u[1] - u[3] / p[5] - u[3] / p[6]
end

## ODEProblem definition with temporary u0 and parameters
u0 = [0.1, 0, 0]
prob = ODEProblem(rate_equations, u0, tspan, lb)

## definition loss function
function loss(par)
    l = 0.0
    function prob_func(prob, i, repeat)
        remake(prob, u0 = [par[num_pars+i], 0, 0])
    end
    ensemble_simu = EnsembleProblem(prob, prob_func = prob_func)
    sol = Array(solve(ensemble_simu, Rosenbrock23(), p = par[1:num_pars-1], EnsembleThreads(),
		              trajectories = num_traj, saveat = t,save_idxs = [1],
		              reltol=1e-6,abstol=1e-8))[1,:,:]
    for i = 1:num_traj
        model = ((C_scaling * par[num_pars] / ND[i]) .* sol[:,i] .* (sol[:,i] .+ 1)) .+ background[i]
        idx_used = idx_IPL[:,i] .& (model .> 0) # we ignore zeros for the log
        l += sum((model[idx_used]) .- (IPL[idx_used,i] .* log.(model[idx_used])))
    end
    return l
end

## definition of loss function for NLopt
function loss_nlopt(par::Vector,grad::Vector)
    if length(grad) > 0
	    grad[:] = ForwardDiff.gradient(loss,par)
    end
    return loss(par)
end

## multi run optimization
function multirun_opt(n_iter,n_runs)
    results = Array{Any}(undef,(n_runs+2,dims+3))
    results[1:n_runs,1] = ["Run $i" for i = 1:n_runs]
    results[end-1:end,[1,2,dims+3]] = ["LB" "x" "x";"UP" "x" "x"]
    results[end-1:end,3:end-1] = hcat(lb,ub)'
    for i = 1:n_runs
        seed_num = abs(rand(Int))
        println("$(seed_num) RNG Seed")
        Random.seed!(seed_num)
        println("Run $(i)")
        flush(stdout)
        optim_setup = bbsetup(loss,Method = :adaptive_de_rand_1_bin_radiuslimited,
                              SearchRange = collect(zip(lb, ub)),
                              NumDimensions = dims, MaxFuncEvals = n_iter,
                              TraceMode = :compact, TraceInterval = 30)
        optim_result = bboptimize(optim_setup)
        p_optim = best_candidate(optim_result)
        optl = Opt(:LD_LBFGS,dims)
        optl.min_objective = loss_nlopt
        optl.ftol_abs = ftol_loc_opt
        optl.maxtime = max_time_loc_opt
        (optfl,optpl,retl) = optimize(optl,p_optim)
        display(optfl)
        println("")
        display(optpl)
        println("")
        results[i,2:end] = vcat(optfl,optpl,Int(seed_num))
    end
    return results
end
## multi run optimization
res_opt = multirun_opt(n_iter,n_runs)

## output preparation optimization results
df_opt = DataFrame(res_opt,:auto)
var_names_opt = vcat(["Run/Bounds","F_min", "τ_r0", "τ_nr0", "τ_2t",
		              "τ_1t", "τ_1e", "τ_1d", "α", "β", "r","C"],
		              [repeat("η_0$i", 1) for i = 1:num_traj],"Seed")
rename!(df_opt, var_names_opt)

## output results
CSV.write("$(name_measurement)_opt_par.csv", df_opt)

## message
println("Script has finished")
## path to TRPL folder --> set it before running the code!
path_to_trpl = "D:\\Projects\\TRPL" # some path to the TRPL folder (just an example here)

## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")

## loading packages
using OrdinaryDiffEq, LinearAlgebra, BlackBoxOptim, ForwardDiff,
	  Printf, Distributions, DataFrames, CSV, XLSX, NLopt

## directory
cd(joinpath(path_to_trpl,"Data\\Real"))

## output results? (true or false)
out = false

## parameters measurement data
name_measurement = "sample1"
N_DA = 5e16 # sample1: N_DA = 5e16, sample2: N_DA = 1e17

## load measurement data
data_measurement = Matrix(CSV.read(name_measurement*".csv",DataFrame,skipto=2))

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

## definition rate equations
function rate_equations(du, u, p, t)
    du[1] = -((1 - u[2]) / p[3]) * u[1] - ((1 - u[3]) / p[4]) * u[1] + (1 / p[5]) * u[3] * 1 / p[7] -
            ((1 + u[1]) / p[1]) * u[1] - ((1 + u[1]) / (p[2] * (1 + p[9] * u[1]))) * u[1]
    du[2] = ((1 - u[2]) / (p[3])) * p[8] * u[1]
    du[3] = ((1 - u[3]) / p[4]) * p[7] * u[1] - u[3] / p[5] - u[3] / p[6]
end

## some inital parameters
# 		         τ_r0,  τ_nr0, τ_d, τ_0,  τ_1, τ_l,  N_l,  N_D,    C,  r , η_0s
init_par = vcat([   1,      1,   1,   1,    1,   1,  1e8,  1e8,  1e4, 10], repeat([1e-6],num_traj))
## ODEProblem definition with temporary u0 and parameters
u0 = [0.1, 0, 0]
prob = ODEProblem(rate_equations, u0, tspan, init_par)

## definition loss function
function loss(par)
    l = 0.0
    function prob_func(prob, i, repeat)
        remake(prob, u0 = [par[num_pars+i], 0, 0])
    end
    ensemble_simu = EnsembleProblem(prob, prob_func = prob_func)
    sol = Array(solve(ensemble_simu, Rosenbrock23(), p = par[1:num_pars-1], EnsembleThreads(),
		      		  trajectories = num_traj, saveat = t,save_idxs = [1],
		      		  reltol = 1e-6,abstol = 1e-8))[1,:,:]
    for i = 1:num_traj
		model = ((C_scaling * par[num_pars] / ND[i]) .* sol[:,i] .* (sol[:,i] .+ 1)) .+ background[i]
		idx_used = idx_IPL[:,i] .& (model .> 0)
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

## load all optimization results
cd(joinpath(path_to_trpl,"Optimizations\\Real\\",name_measurement))
opt_outputs = Vector{String}()
for (root, dirs, files) in walkdir(pwd())
    for file in files
		if file[end-3:end] == ".csv"
			push!(opt_outputs,joinpath(root,file))
		end
    end
end

## process optimization parameters and results (needed for output)
n_runs = size(Matrix(CSV.read(opt_outputs[1],DataFrame,skipto=2)),1)-2
n_jobs = length(opt_outputs)
n_opt = Int(n_jobs * n_runs)
println("$n_opt total optimizations: $n_jobs batchjobs with $n_runs optimizations each")
param = Array{Any}(undef,(n_opt,dims))
mle_loss = Array{Any}(undef,(n_opt,1))
rand_seeds = Array{Any}(undef,(n_opt,1))
for i = 0:n_jobs-1
	param[1+(n_runs*i):n_runs+(n_runs*i),:] = Matrix(CSV.read(opt_outputs[i+1],DataFrame,skipto=2))[1:n_runs,3:end-1]
	mle_loss[1+(n_runs*i):n_runs+(n_runs*i)] = Matrix(CSV.read(opt_outputs[i+1],DataFrame,skipto=2))[1:n_runs,2]
	rand_seeds[1+(n_runs*i):n_runs+(n_runs*i)] = Matrix(CSV.read(opt_outputs[i+1],DataFrame,skipto=2))[1:n_runs,end]
end
lb = Matrix(CSV.read(opt_outputs[1],DataFrame,skipto=2))[end-1,3:end-1]
ub = Matrix(CSV.read(opt_outputs[1],DataFrame,skipto=2))[end,3:end-1]
rand_seeds = Int.([parse(BigFloat,rand_seeds[i]) for i = 1:n_opt])
mle_loss = [parse(Float64,mle_loss[i]) for i = 1:n_opt]
opt_results = hcat(mle_loss,param)
sorted_opt_results = opt_results[sortperm(opt_results[:, 1]),:]

## get best fit parameters
p_fit = sorted_opt_results[1,2:end]

## check that the loss function defined here gives the same result
println("sanity check for loss: loss(p_fit) - loss(p_csv) = "*
        "$(loss(p_fit) .- sorted_opt_results[1,1])")

## rescaling parameters for differentiation
p_r = copy(p_fit)
sf = 10 .^ -floor.(log10.(p_r)) |> Diagonal
p_rescaled = sf * p_r 

## define rescaled loss function
function loss_rescaled(par)
	_par = inv(sf) * par
    l = 0.0
    function prob_func(prob, i, repeat)
        remake(prob, u0 = [_par[num_pars+i], 0, 0])
    end
    ensemble_simu = EnsembleProblem(prob, prob_func = prob_func)
    sol = Array(solve(ensemble_simu, Rosenbrock23(), p = _par[1:num_pars-1], EnsembleThreads(),
		     		  trajectories = num_traj, saveat = t,save_idxs = [1],
		      		  reltol = 1e-6,abstol = 1e-8))[1,:,:]
    for i = 1:num_traj
		model = ((C_scaling * _par[num_pars] / ND[i]) .* sol[:,i] .* (sol[:,i] .+ 1)) .+ background[i]
		idx_used = idx_IPL[:,i] .& (model .> 0)
		l += sum((model[idx_used]) .- (IPL[idx_used,i] .* log.(model[idx_used])))
    end
    return l
end

## define rescaled loss for NLopt
function loss_rescaled_nlopt(par::Vector,grad::Vector)
    if length(grad) > 0
        grad[:] = ForwardDiff.gradient(loss_rescaled,par)
    end
    display(loss_rescaled(par))
    return loss_rescaled(par)
end

## sanity check for rescaled loss function
println("sanity check for loss: loss_rescaled(p_rescaled) - loss(p_csv) = "*
        "$(loss_rescaled(p_rescaled) .- sorted_opt_results[1,1])")

## local optimization without bounds around the found optimum (LBFGS: local optimization)
opt = Opt(:LD_LBFGS,dims)
opt.min_objective = loss_rescaled_nlopt
opt.ftol_abs = 1e-14
opt.maxtime = 600
(optf,p_opt,ret) = optimize(opt,p_rescaled)

## calculate gradient
g = ForwardDiff.gradient(loss_rescaled,p_opt)

## calculate hessian
H = ForwardDiff.hessian(loss_rescaled,p_opt)

## calculate eigenvalues of the hessian
e = eigvals(H)

## calculate standard errors of fit parameters
se = inv(sf) * inv(H) * inv(sf) |> diag .|> sqrt

## convert parameters N_st, and N_dt back to ffective DOS 
p_fit_output = copy(p_fit) 
p_fit_output[7:8] = N_DA ./ p_fit_output[7:8]
p_opt_output = inv(sf) * copy(p_opt) 
p_opt_output[7:8] = N_DA ./ p_opt_output[7:8]
se_output = copy(se)
se_output[7:8] = se_output[7:8] .* p_opt_output[7:8] .^2 ./ N_DA

## parameter table
par_table = Array{Any}(undef,7,dims+2)
par_table[1,1] = "p_fit"
par_table[1,2:end] = vcat(loss(p_fit),p_fit_output)
par_table[2,1] = "p_opt"
par_table[2,2:end] = vcat(loss(inv(sf)*p_opt),p_opt_output)
par_table[3,1] = "ci95"
par_table[3,2:end] = vcat(["/"],2.98*se_output)
par_table[4,1] = "lb"
par_table[4,2:end] = vcat(["/"],lb)
par_table[5,1] = "ub"
par_table[5,2:end] = vcat(["/"],ub)
par_table[6,1] = "gradient"
par_table[6,2:end] = vcat(["/"],g)
par_table[7,1] = "hessian eigenvalues (sorted)"
par_table[7,2:end] = vcat(["/"],e)

## optimation result table
opt_table = Array{Any}(undef,n_opt,dims+3)
opt_table[:,2:end] = hcat(opt_results,rand_seeds)
opt_table[:,1] = ["Run $i" for i = 1:n_opt]

## data for output
data_table = Array{Any}(undef,size(data_measurement,1)+1,num_traj+1)
data_table[1,1] = "doping"
data_table[1,2:end] = vcat(N_DA,repeat(["/"],num_traj-1))
data_table[2:end,:] = data_measurement
data_table[2,1] = "power"
data_table[3,1] = "ND"
data_table[4,1] = "background"

## variable names 
var_names = vcat(["","l_opt", "τ_r0", "τ_nr0", "τ_d",
				  "τ_0", "τ_1", "τ_l", "N_l", "N_D", "r", "C"],
	 			  [repeat("η_0$i", 1) for i = 1:num_traj])

## generate output data
df_res = DataFrame(par_table,:auto)
rename!(df_res,var_names)
df_res 
df_opt = DataFrame(opt_table,:auto)
rename!(df_opt,vcat(var_names,"Seed"))
df_data = DataFrame(data_table,:auto)
rename!(df_data,vcat("/",["Transient $i" for i = 1:num_traj]))

## output
if out
	## output file path
	file_name = "$name_measurement.xlsx"
	xlsx_path = joinpath(path_to_trpl,"Results\\Real",file_name)

	## save output
	if isfile(xlsx_path)
		rm(xlsx_path)
	end
	XLSX.writetable(xlsx_path, Results = (collect(DataFrames.eachcol(df_res)), DataFrames.names(df_res)),
					Opts = (collect(DataFrames.eachcol(df_opt)), DataFrames.names(df_opt)),
					Data = (collect(DataFrames.eachcol(df_data)), DataFrames.names(df_data)))
end
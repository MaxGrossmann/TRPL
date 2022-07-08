## path to TRPL folder --> set it before running the code!
path_to_trpl = "D:\\Projects\\TRPL" # some path to the TRPL folder (just an example here)

## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")

## loading packages
using OrdinaryDiffEq, LinearAlgebra, BlackBoxOptim, ForwardDiff, DelimitedFiles,
	  Printf, Distributions, DataFrames, CSV, XLSX, NLopt

## directory
cd(joinpath(path_to_trpl,"Data\\Syn\\CSV"))

## output results? (true or false)
out = false

## index of synthetic data set
idx_syn_data = 1

## confidence level (1-α)
alpha = 0.05 # corresponds to a 95% confidence level
q = quantile(Normal(),1-alpha/2)


## load measurement data
data_measurement = Matrix(CSV.read("syn_data_$idx_syn_data.csv",DataFrame,skipto=2))

## function for data processing (data structre is explained in README.md)
function data_processing(data)
	t = data[4:end,1]
	num_traj = size(data,2)-1
	num_pars = 10
	dims = num_pars + num_traj
	C_scaling = data[1,2]
	ND = 10 .^ data[2,2:end]
	background = data[3,2:end]
	IPL = data[4:end,2:num_traj+1]
	tspan = (t[1], t[end])
	N = length(t)
	idx_IPL = zeros(Bool,size(IPL))
	for i = 1:num_traj
		idx_IPL[:,i] = IPL[:,i] .> 0
	end
    return t, tspan, num_traj, num_pars, dims, IPL, idx_IPL, ND, background, N, C_scaling
end

## process synthetic data
t, tspan, num_traj, num_pars, dims, IPL, idx_IPL, ND, background, N, C_scaling = data_processing(data_measurement)

## definition rate equations
function rate_equations(du, u, p, t)
    du[1] = -((1 - u[2]) / p[3]) * u[1] - ((1 - u[3]) / p[4]) * u[1] + (1 / p[5]) * u[3] * 1 / p[7] -
            ((1 + u[1]) / p[1]) * u[1] - ((1 + u[1]) / (p[2] * (1 + p[9] * u[1]))) * u[1]
    du[2] = ((1 - u[2]) / (p[3])) * p[8] * u[1]
    du[3] = ((1 - u[3]) / p[4]) * p[7] * u[1] - u[3] / p[5] - u[3] / p[6]
end

## some inital parameters
# 		         τ_r0,  τ_nr0, τ_d, τ_0,  τ_1, τ_l,  N_l,  N_D,    C,  r ,  η_0s
init_par = vcat([   1,      1,   1,   1,    1,   1,  1e8,  1e8,  1e4, 10], repeat([1e-6],num_traj))

## ODEProblem definition with temporary u0 and parameters
u0 = [0.1, 0, 0]
prob = ODEProblem(rate_equations, u0, tspan, init_par)

## function to generate the model for plotting
function generateipl(par)
    function prob_func(prob, i, repeat)
        remake(prob, u0 = [par[num_pars+i], 0, 0])
    end
    ensemble_simu = EnsembleProblem(prob, prob_func = prob_func)
    sol = Array(solve(ensemble_simu, Rosenbrock23(), p = par[1:num_pars-1], EnsembleThreads(),
		      		  trajectories = num_traj, saveat = t, save_idxs = [1],
		      		  reltol = 1e-6, abstol = 1e-8))[1,:,:]
    ipl = zeros(size(sol))
    for i = 1:num_traj
    	ipl[:,i] = ((C_scaling * par[num_pars] / ND[i]) .* sol[:,i] .* (sol[:,i] .+ 1)) .+ background[i]
    end
    return ipl
end

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

## load all optimization results
cd(joinpath(path_to_trpl,"Optimizations\\Syn\\syn_data_$idx_syn_data"))
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
lb = Float64.(Matrix(CSV.read(opt_outputs[1],DataFrame,skipto=2))[end-1,3:end-1])
ub = Float64.(Matrix(CSV.read(opt_outputs[1],DataFrame,skipto=2))[end,3:end-1])
rand_seeds = Int.([parse(BigFloat,rand_seeds[i]) for i = 1:n_opt])
mle_loss = [parse(Float64,mle_loss[i]) for i = 1:n_opt]
opt_results = hcat(mle_loss,param)
sorted_opt_results = opt_results[sortperm(opt_results[:, 1]),:]

## get best fit parameters
p_fit = Float64.(sorted_opt_results[1,2:end])

## check that the loss function defined here gives the same result
println("sanity check for loss: loss(p_fit) - loss(p_csv) = "*
		"$(loss(p_fit) .- sorted_opt_results[1,1])")

## load true parameters
cd(joinpath(path_to_trpl,"Data\\Syn\\Params"))
input_params = readdlm("syn_data_$idx_syn_data.txt",'\n')
p_true = let expr = Meta.parse(input_params[3])
@assert expr.head == :vect
Array(expr.args)
end

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
		      		  reltol=1e-6,abstol=1e-8))[1,:,:]
    for i = 1:num_traj
		model = ((C_scaling * _par[num_pars] / ND[i]) .* sol[:,i] .* (sol[:,i] .+ 1)) .+ background[i]
		idx_used = idx_IPL[:,i] .& (model .> 0) # we ignore zeros for the log
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

## rescaled loss for NLopt
function loss_rescaled_nlopt(par::Vector,grad::Vector)
    if length(grad) > 0
        grad[:] = ForwardDiff.gradient(loss_rescaled,par)
    end
    display(loss_rescaled(par))
    return loss_rescaled(par)
end

## local optimization without bounds around the found optimum (LBFGS: local optimization)
opt = Opt(:LD_LBFGS,dims)
opt.min_objective = loss_rescaled_nlopt
opt.ftol_abs = 1e-14
opt.maxtime = 600
(optf,p_opt,ret) = optimize(opt,p_rescaled)

## calculate gradient
gL = ForwardDiff.gradient(loss_rescaled,p_opt)

## calculate hessian
HL = ForwardDiff.hessian(loss_rescaled,p_opt)

## calculate eigenvalues of the hessian
eL = eigvals(HL)

## function to calculate the rescaled model gradient for the expection value of the hessian matrix (fisher information)
function rescaled_lambda(par,eta0_idx,t_idx)
    _par = inv(sf) * par
    prob = ODEProblem(rate_equations, [_par[num_pars + eta0_idx], 0, 0], tspan, _par[1:num_pars])
    sol = Array(solve(prob, Rosenbrock23(), saveat = t, save_idxs = [1],
		      		  reltol = 1e-6, abstol = 1e-8))[1,:]
    ipl = ((C_scaling * _par[num_pars] / ND[eta0_idx]) .* sol[t_idx] .* (sol[t_idx] .+ 1)) .+ background[eta0_idx]
    return ipl
end

## calculate expection value of the hessian matrix (fisher information)
grad_lambda(p,i,j) = ForwardDiff.gradient(x->rescaled_lambda(x,i,j),p)
H = zeros(dims,dims)
for i = 1:num_traj
    for j = 1:length(t)
        global H += (grad_lambda(p_opt,i,j) * grad_lambda(p_opt,i,j)')/rescaled_lambda(p_opt, i, j)
    end
end

## calculate standard errors of fit parameters
se = inv(Diagonal(sf)) * inv(H) * inv(Diagonal(sf)) |> diag .|> sqrt

## generate fit data for plot
IPL_opt = generateipl(inv(sf) * p_opt)

## prepare output
p_true_output = copy(p_true) 
p_fit_output = copy(p_fit) 
p_opt_output = inv(sf) * copy(p_opt) 
se_output = copy(se)

## error table
# p_fit is best global optimization result
# p_opt is the final result after the reparameterized local optimization (started at p_fit)
par_table = Array{Any}(undef,11,dims+2)
par_table[1,1] = "p_fit"
par_table[1,2:end] = vcat(loss(p_fit),p_fit_output)
par_table[2,1] = "p_opt"
par_table[2,2:end] = vcat(loss(inv(sf)*p_opt),p_opt_output)
par_table[3,1] = "p_true"
par_table[3,2:end] = vcat(loss(p_true),p_true_output)
par_table[4,1] = "difference"
par_table[4,2:end] = vcat(loss(inv(sf)*p_opt) .- loss(p_true),abs.(p_opt_output .- p_true_output))
par_table[5,1] = "rel. error in %"
par_table[5,2:end] = vcat(abs.(1 .- (loss(inv(sf)*p_opt) ./ loss(p_true))) * 100,
						    abs.(1 .- (inv(sf)*p_opt ./ p_true)) * 100)
par_table[6,1] = "avg. rel. error in %"
par_table[6,2:end] = vcat(["/"],mean(abs.(1 .- (inv(sf)*p_opt ./ p_true)) * 100),
				            repeat(["/"],dims-1))
par_table[7,1] = "cihw95"
par_table[7,2:end] = vcat(["/"],q*se_output)
par_table[8,1] = "lb"
par_table[8,2:end] = vcat(["/"],lb)
par_table[9,1] = "ub"
par_table[9,2:end] = vcat(["/"],ub)
par_table[10,1] = "gradient"
par_table[10,2:end] = vcat(["/"],gL)
par_table[11,1] = "hessian eigenvalues (sorted)"
par_table[11,2:end] = vcat(["/"],eL)

## optimation result table
opt_table = Array{Any}(undef,n_opt,dims+3)
opt_table[:,2:end] = hcat(opt_results,rand_seeds)
opt_table[:,1] = ["Run $i" for i = 1:n_opt]

## data for output
data_table = Array{Any}(undef,size(data_measurement,1),num_traj+1)
data_table[:,:] = data_measurement
data_table[1,1] = "C_scaling"
data_table[1,3:end] = repeat(["/"],num_traj-1)
data_table[2,1] = "ND"
data_table[3,1] = "background"

## variable names 
var_names = vcat(["","l_opt", "τ_r0", "τ_nr0", "τ_2t",
				  "τ_1t", "τ_1e", "τ_1d", "α", "β", "r", "C"],
	 			  [repeat("η_0$i", 1) for i = 1:num_traj])

## generate output data
df_res = DataFrame(par_table,:auto)
rename!(df_res,var_names)
df_fit = DataFrame(hcat(t,IPL_opt),:auto)
rename!(df_fit,vcat("t",["Transient $i" for i = 1:num_traj]))
df_opt = DataFrame(opt_table,:auto)
rename!(df_opt,vcat(var_names,"Seed"))
df_data = DataFrame(data_table,:auto)
rename!(df_data,vcat("/",["Transient $i" for i = 1:num_traj]))

## outputs
if out 
	## output file path
	file_name = ("syn_data_$idx_syn_data.xlsx")
	xlsx_path = joinpath(path_to_trpl,"Results\\Syn",file_name)
	
	## save output
	if isfile(xlsx_path)
		rm(xlsx_path)
	end
	XLSX.writetable(xlsx_path, Results = (collect(DataFrames.eachcol(df_res)), DataFrames.names(df_res)),
					Fits = (collect(DataFrames.eachcol(df_fit)), DataFrames.names(df_fit)),
					Opts = (collect(DataFrames.eachcol(df_opt)), DataFrames.names(df_opt)),
					Data = (collect(DataFrames.eachcol(df_data)), DataFrames.names(df_data)))
end

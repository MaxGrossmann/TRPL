## path to TRPL folder --> set it before running the code!
path_to_trpl = "C:\\TRPL" # some path to the TRPL folder
## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")
## loading packages
using OrdinaryDiffEq, LinearAlgebra, StaticArrays, Random, NLopt, DelimitedFiles,
	  BlackBoxOptim, Printf, PoissonRandom, DataFrames, Plots, CSV
## backend quick plots
pyplot()
default(size = (800, 600))
## save?
save_output = true
## noise? 
noise = true
## experimental parameters 
background = [1, 2, 12, 26, 18]
ND         = [1, 1,  2,  3,  4]
C_scaling  = 1e6
## parameters
# 		       τ_r0, τ_nr0,  τ_d, τ_0,    τ_1,  τ_l,     α,     β, r,      C 
parameters = [100.0,  75.0, 15.0, 5.0, 1300.0, 85.0, 100.0, 105.0, 5, 1250.0]
eta0s      = [5.0e-5, 0.0005, 0.005, 0.5, 0.5]
## time
tspan = [0.0, 500.0]
dt    = 0.25
## number of data set
idx_syn = 4
## useful variables
num_pars = 10
num_traj = length(eta0s)
dims =  num_pars + num_traj
## output name
name_syn_data = "syn_data_$(idx_syn)"
## definition rate equations
function rate_equations(du, u, p, t)
    du[1] = -((1 - u[2]) / p[3]) * u[1] - ((1 - u[3]) / p[4]) * u[1] + (1 / p[5]) * u[3] * 1 / p[7] -
            ((1 + u[1]) / p[1]) * u[1] - ((1 + u[1]) / (p[2] * (1 + p[9] * u[1]))) * u[1]
    du[2] = ((1 - u[2]) / (p[3])) * p[8] * u[1]
    du[3] = ((1 - u[3]) / p[4]) * p[7] * u[1] - u[3] / p[5] - u[3] / p[6]
end
## ODEProblem definition with temporary u0 and parameters
u0 = [0.1, 0, 0]
t = [tspan[1]:dt:tspan[2];]
p_true = vcat(parameters,eta0s)
prob = ODEProblem(rate_equations, u0, tspan, p_true)
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
    	ipl[:,i] = ((C_scaling * par[num_pars] / (10 ^ ND[i])) .* sol[:,i] .* (sol[:,i] .+ 1)) .+ background[i]
    end
    return ipl
end
## get synthetic data with poisson noise if wanted
IPL = generateipl(p_true)
if noise
	rng_seed = abs(rand(Int))
	IPL = pois_rand.(Random.seed!(rng_seed),IPL)
end
## calculate loss goal
idx_IPL = zeros(Bool,size(IPL))
for i = 1:num_traj
	idx_IPL[:,i] = IPL[:,i] .> 0
end
function rate_equations(du, u, p, t)
	du[1] = -((1 - u[2]) / p[3]) * u[1] - ((1 - u[3]) / p[4]) * u[1] + (1 / p[5]) * u[3] * 1 / p[7] -
			((1 + u[1]) / p[1]) * u[1] - ((1 + u[1]) / (p[2] * (1 + p[9] * u[1]))) * u[1]
	du[2] = ((1 - u[2]) / (p[3])) * p[8] * u[1]
	du[3] = ((1 - u[3]) / p[4]) * p[7] * u[1] - u[3] / p[5] - u[3] / p[6]
end
u0 = [0.1, 0, 0]
init_par = vcat([   1,      1,   1,   1,    1,   1,  1e8,  1e8,  1e4, 10], repeat([1e-6],num_traj))
prob = ODEProblem(rate_equations, u0, tspan, init_par)
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
		model = ((C_scaling * par[num_pars] / (10 ^ ND[i])) .* sol[:,i] .* (sol[:,i] .+ 1)) .+ background[i]
		idx_used = idx_IPL[:,i] .& (model .> 0)
		l += sum((model[idx_used]) .- (IPL[idx_used,i] .* log.(model[idx_used])))
	end
	return l
end
l = loss(p_true)
println("loss goal = $l") 
## plotting font size
fs = 18
## plot data
ylim = (10^0,10^(ceil(log10(maximum(IPL)))))
plot_window = plot(palette = palette(:seaborn_colorblind, num_traj), xlims = (t[1], t[end]), yscale=:log10,
				   ylims= ylim, title = "Synthetic data set $idx_syn", titlefontsize = fs, grid = true,
				   xguidefontsize = fs, yguidefontsize = fs, dpi = 600, framestyle = :box, xtickfontsize = fs,
				   ytickfontsize = fs, xlabel = "t/ns", ylabel = "Photon Counts", legendfontsize = fs-4, leg = :topright,
				   margin = 10Plots.mm)
for i = num_traj:-1:1
	plot!(t, IPL[:,i],linealpha=0.75, linewidth=1.5, label=@sprintf("η₀ = %.2E", p_true[num_pars+i]))
end
display(plot_window)
if save_output
	## save input parameters
	io = open(string(path_to_trpl,"\\Data\\Syn\\Params\\$name_syn_data.txt"), "w")
		println(io, "parameters")
		println(io, "τ_r0,  τ_nr0,  τ_d,  τ_0,   τ_1,  τ_l,   α,  β,  r,    C, η0s")
		println(io, p_true)
		println(io, "C scaling")
		println(io, C_scaling)
		println(io, "background")
		println(io, background)
		println(io, "ND")
		println(io, ND)
		println(io, "time span")
		println(io, tspan)
		println(io, "dt")
		println(io, dt)
		println(io, "loss goal")
		println(io, l)
		if noise
			println(io, "rng seed")
			println(io, rng_seed)
		end
	close(io)
	## save data for fit
	output = Array{Any}(undef,(size(IPL,1)+3,num_traj+1))
	output[1,2] = C_scaling
	output[2,2:end] = ND
	output[3,2:end] = background
	output[4:end,1] = t
	output[4:end,2:end] = IPL
	output[1:3,1] .= 0
	output[1,3:end] .= 0
	df = DataFrame(output,:auto)
	CSV.write(string(path_to_trpl,"\\Data\\Syn\\CSV\\$name_syn_data.csv"), df)
end
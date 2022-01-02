## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")
## loading packages
using CSV, DataFrames, XLSX, Printf, PrettyTables
## parameters measurement data
idx_syn_data = 5
## load measurement data
cd(string(homedir(),"\\OneDrive\\Uni_Master\\TRPL\\Data\\Syn\\data"))
data_measurement = Matrix(CSV.read("syn_data_$idx_syn_data.csv",DataFrame,skipto=2))
## function for data processing
function data_processing(data)
	t = data[4:end,1]
	num_traj = size(data,2)-1
	num_pars = 10
	dims = num_pars + num_traj
	N_DA = data[1,2]
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
    return t, tspan, num_traj, num_pars, dims, IPL, idx_IPL, N_DA, ND, background, N, C_scaling
end
## loading of experimental data
t, tspan, num_traj, num_pars, dims, IPL, idx_IPL, N_DA, ND, background, N, C_scaling = data_processing(data_measurement)
## directory
cd(string(homedir(),"\\OneDrive\\Uni_Master\\TRPL\\Data\\Final\\Supplement\\data"))
## get xlsx and sort them
function custom_cmp(x::String)
    number_idx = findfirst(isdigit, x)
    str, num = SubString(x, 1, number_idx-1), SubString(x, number_idx, length(x)-5)
    return str, parse(Int, num)
 end
files =  sort(filter(x->endswith(x, ".xlsx"), readdir()), by = custom_cmp)
## get results sheet
println(files[idx_syn_data])
println("")
xf = XLSX.readxlsx(files[idx_syn_data])
sh = xf["Results"]
num_pars = 10
## get parameters
data = Array{Any}(undef,5,num_pars + 1)
data[1,1] = ""
data[1,2:end] = ["\$\\phi \\tau_r^0\$", "\$\\tau_{nr}^0\$" ,"\$\\tau_{dt}\$", "\$\\tau_{st}\$", "\$\\tau_{se}\$", "\$\\tau_{sd}\$",
                 "\$\\alpha \$",  "\$\\beta\$", "\$r\$", "\$C\$"]
data[2:end,1] = ["\$p_{\\text{true}}\$", "\$p_{\\text{fit}}\$", "\$\\text{SE}\$", "\$ \\delta \$ (%)"]
data[2,2:end] = sh["C4:L4"]
data[3,2:end] = sh["C3:L3"]
data[4,2:end] = sh["C8:L8"]
data[5,2:end] = sh["C6:L6"]
data[2:4,end] .*= C_scaling
## format numbers
sf = floor.(Int,log10.(data[2:end,2:end]))
data_r = round.(data[2:end,2:end], sigdigits = 3)
data[2:end-1,2:9] = [@sprintf("\$ %.2f\$",data_r[i,j]) for i = 1:3, j = 1:8]
data[2:end-1,end-1] = [@sprintf("\$ %.2f \$", data[i,end-1]) for i = 2:4]
data[2:end-1,end] = [@sprintf("\$ %.2f \\cdot 10^{%.i} \$", data_r[i,end] / (10.0 ^ sf[i,end]), sf[i,end]) for i = 1:3]
data[end,2:end] = [@sprintf("\$ %.2f \$",data_r[end,j]) for j = 1:num_pars]
## reorder columns 
table = @view data[:,[1,2,3,10,4,5,6,7,8,9,11]]
## tex table
pretty_table(table, backend = Val(:latex), table_type = :tabular, noheader = true,
             alignment=:c, body_hlines = [1],  hlines = :none)
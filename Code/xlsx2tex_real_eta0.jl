## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")
## loading packages
using XLSX, Printf, PrettyTables
## directory
cd(string(homedir(),"\\OneDrive\\Uni_Master\\TRPL\\Data\\Final\\Real\\data"))
## index of synthetic data set
idx = 1
## get xlsx and sort them
files =  filter(x->endswith(x, ".xlsx"), readdir())
## get results sheet
println(files[idx])
println("")
xf = XLSX.readxlsx(files[idx])
sh = xf["Results"]
num_traj = size(sh[:],2) - 2 - 10
## get η₀
data = Array{Any}(undef,3,num_traj + 1)
data[1,1] = ""
data[1,2:end] = ["\$\eta_0^$i\$" for i = 1:5]
data[2:end,1] = ["\$\eta_{\text{fit}}\$", "\$\text{SE}\$"]
data[2:3,2:end] = sh["M3:Q4"]
## format numbers
sf = floor.(Int,log10.(data[2:end,2:end]))
data_r = round.(data[2:end,2:end] ./ 10.0 .^ sf, sigdigits = 3)
data[2:end,2:end] = [@sprintf("\$ %.2f \\cdot 10^{%i}\$",data_r[i,j],sf[i,j]) for i = 1:2, j = 1:num_traj]
## tex table
pretty_table(data, backend = Val(:latex), table_type = :tabular, noheader = true,
             alignment=:c, body_hlines = [1],  hlines = :none)
## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")
## loading packages
using XLSX, Printf, PrettyTables
## directory
cd(string(homedir(),"\\OneDrive\\Uni_Master\\TRPL\\Data\\Final\\Supplement\\data"))
## index of synthetic data set
idx = 5
## get xlsx and sort them
function custom_cmp(x::String)
    number_idx = findfirst(isdigit, x)
    str, num = SubString(x, 1, number_idx-1), SubString(x, number_idx, length(x)-5)
    return str, parse(Int, num)
 end
files =  sort(filter(x->endswith(x, ".xlsx"), readdir()), by = custom_cmp)
## get results sheet
println(files[idx])
println("")
xf = XLSX.readxlsx(files[idx])
sh = xf["Results"]
num_traj = size(sh[:],2) - 2 - 10
## get η₀
data = Array{Any}(undef,5,num_traj + 1)
data[1,1] = ""
data[1,2:end] = ["\$\eta_0^$i\$" for i = 1:num_traj]
data[2:end,1] = ["\$\eta_{\text{true}}\$", "\$\eta_{\text{fit}}\$", "\$\text{SE}\$", "\$ \\delta \$ (%)"]
cols = ["N", "O", "P", "Q", "R", "S", "T", "U", "V"]
data[2,2:end] = sh["M4:$(cols[num_traj-1])4"]
data[3,2:end] = sh["M3:$(cols[num_traj-1])3"]
data[4,2:end] = sh["M8:$(cols[num_traj-1])8"]
data[5,2:end] = sh["M6:$(cols[num_traj-1])6"]
## format numbers
sf = floor.(Int,log10.(data[2:end,2:end]))
data_r = round.(data[2:end,2:end] ./ 10.0 .^ sf, sigdigits = 3)
data[2:end-1,2:end] = [@sprintf("\$ %.2f \\cdot 10^{%i}\$",data_r[i,j],sf[i,j]) for i = 1:3, j = 1:num_traj]
data[end,2:end] = [@sprintf("\$ %.2f \$",data_r[end,j] .* 10.0 .^ sf[end,j]) for j = 1:num_traj]
## tex table
pretty_table(data, backend = Val(:latex), table_type = :tabular, noheader = true,
             alignment=:c, body_hlines = [1],  hlines = :none)
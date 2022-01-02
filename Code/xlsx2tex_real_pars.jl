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
num_pars = 10
## get parameters
data = Array{Any}(undef,3,num_pars + 1)
data[1,1] = ""
data[1,2:end] = ["\$\\phi \\tau_r^0\$", "\$\\tau_{nr}^0\$" ,"\$\\tau_{dt}\$", "\$\\tau_{st}\$", "\$\\tau_{se}\$", "\$\\tau_{sd}\$",
                   "\$N_{st}\$",  "\$N_{dt}\$", "\$r\$" ,"\$C\$" ]
data[2:end,1] = ["\$p_{\\text{fit}}\$", "\$\\text{SE}\$"]
data[2:3,2:end] = sh["C3:L4"]
## format numbers
sf = floor.(Int,log10.(data[2:end,2:end]))
data_r = round.(data[2:end,2:end], sigdigits = 3)
data[2:end,2:end-4] = [@sprintf("\$ %.2f\$",data_r[i,j]) for i = 1:2, j = 1:6]
data[2:end,end-3:end-2] = [@sprintf("\$ %.2f \\cdot 10^{%.i} \$", data_r[i,j] / (10.0 ^ sf[i,j]), sf[i,j]) for i = 1:2, j = 7:8]
data[2:end,end-1] = [@sprintf("\$ %.2f \$", data[i,end-1]) for i = 2:3]
data[2:end,end] = [@sprintf("\$ %.2f \\cdot 10^{%.i} \$", data_r[i,end] / (10.0 ^ sf[i,end]), sf[i,end]) for i = 1:2]
## reorder columns 
table = @view data[:,[1,2,3,10,4,5,6,7,8,9,11]]
## tex table
pretty_table(table, backend = Val(:latex), table_type = :tabular, noheader = true,
             alignment=:c, body_hlines = [1],  hlines = :none)
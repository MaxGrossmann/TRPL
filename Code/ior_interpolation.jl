## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")
## percentage x of Al, wavelengths λ and index of refraction n from publication
x = [0.491 0.590; 0.7 0.804]
lambda = [0.405 0.520]
n = [4.5889 + 0.83390im 3.8264 + 0.19880im;
     4.4384 + 0.68930im 3.7477 + 0.15232im;
     4.2759 + 0.53061im 3.6807 + 0.05888im;
     4.1150 + 0.40503im 3.5071 + 0.00385im]
N = real(n)
k = imag(n)
## interpolation of n for wanted x
x_wanted = [0.509 0.720]
l = length(lambda)
m = length(x_wanted)
N_wanted = zeros(m,l)
k_wanted = zeros(m,l)
for i = 1:l
    for j = 0:m-1
        local coeffs_N = [x[j+1,:] ones(2)]\N[(1+j*2):(1+j*2)+1,i]
        local coeffs_k = [x[j+1,:] ones(2)]\k[(1+j*2):(1+j*2)+1,i]
        N_wanted[i,j+1] = coeffs_N[1]*x_wanted[j+1] + coeffs_N[2]
        k_wanted[i,j+1] = coeffs_k[1]*x_wanted[j+1] + coeffs_k[2]
    end
end
N_wanted = round.(N_wanted,sigdigits=4)
k_wanted = round.(k_wanted,sigdigits=4)
## output solution
println("")
for j = 1:m
    for i = 1:l
    println("n(x=$(x_wanted[j]), λ = $(lambda[i])) = $(N_wanted[i,j]) + $(k_wanted[i,j])i")
    end
end

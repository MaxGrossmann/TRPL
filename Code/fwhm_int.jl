## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")

## loading packages
using NumericalIntegration, Plots

## normalized gaussian 
n = 5001
x = Array(range(-100,stop=100,length=n))
y = 1/sqrt(2*pi) * exp.(-x .^ 2 / 2)

## check normalization
println("Sanity Check: full integral = $(integrate(x,y))")

## calculate FWHM positions
fwhm1 = sortperm(abs2.(y.- (maximum(y[1:floor(Int,n/2)])/2)))[1]
fwhm2 = n - fwhm1 + 1 

## plot gaussian and mark FWHMs
plt = plot(x,y,label=nothing,linecolor=:black,xlims=(-20,20))
plot!([x[fwhm1]],[y[fwhm1]],label=nothing,marker=:circle,markercolor=:orange)
plot!([x[fwhm2]],[y[fwhm2]],label=nothing,marker=:circle,markercolor=:orange)
display(plt)

## output results
println("area of FWHM: $(round(integrate(x[fwhm1:fwhm2],y[fwhm1:fwhm2]),sigdigits=3))")
## activate package environment
cd(@__DIR__)
using Pkg
Pkg.activate("Project.toml")

## parameters (calculation for perpendicular incident)
λ = 520 # nm
N = [1.0; 3.812 + 0.1903im; 4.1891+0.3633im] # air - window layer - absorber layer
D = 100 # thickness of window layer
# (refractive indices are from https://refractiveindex.info/)

## function for fresnel coefficients (for perpendicular incident)
function fresnel_coeff(n,m,pol)
    if pol == "s"
        r = (n-m)/(n+m)
        t = (2*n)/(n+m)
    elseif pol == "p"
        r = (m-n)/(m+n)
        t = (2*n)/(m+n)
    end
    return r,t
end

## function for interface propagation matrix
function interface_prop(n,m,pol)
    r,t = fresnel_coeff(n,m,pol)
    I = 1/t * [1 r;r 1]
    return I
end

## function for propagation through film Matrix
function film_prop(d,n)
    beta = 2*pi*d*n/λ
    L = [exp(-1im*beta) 0; 0 exp(1im*beta)]
    return L
end

## calcuation of scattering matrix
# Emanuele Centurioni: Generalized matrix method for calculation of internal light energy flux
# in mixed coherent and incoherent multilayers
Ss = interface_prop(N[1],N[2],"s") * film_prop(D,N[2]) * interface_prop(N[2],N[3],"s")
Sp = interface_prop(N[1],N[2],"p") * film_prop(D,N[2]) * interface_prop(N[2],N[3],"p")

## backside transmittance
Ts = abs2(1/Ss[1,1])*real(N[3])/real(N[1])
Tp = abs2(1/Ss[1,1])*real(conj(N[3]))/real(N[1])
T  = round((Ts+Tp)/2,sigdigits=4)

## ouput result
println("")
println("transmission into absorber T = $T")

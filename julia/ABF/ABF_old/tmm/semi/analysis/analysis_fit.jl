using MAT
using Random
using Plots
using Distributions
gr()
include("xi_abf2.jl")


dir = "/Users/pcs/Data/ABF/Rough/rough1"
cd(dir)
#read all the data and the parameter
weak_file = ("config-weak-xi_all.mat")
all_file = ("config-xi_all.mat")
ξ,params = ξ_read(all_file)
ξ_weak, params_weak = ξ_read(weak_file)
# create paraameter space from the configuration
(E, W, θI, θII) =params_expand(params)

ξ_weak_diag = [ξ_weak[1,1,i,i] for i in 1:params.θI_num]

rng = MersenneTwister(12234)
W = 0.00001
a = W*(rand(rng,100000,6).-0.5);

storage = Array{Float64}(undef, 100)
θ = LinRange(0,π/2, 100)
for (i,θi) in enumerate(θ)
    e1 = a[:,1]sin(θi)^2 + a[:,2]cos(θi)^2
    #e2 = a[:,3] .* a[:,4] ./ (a[:,3]sin(θi)^2 .+ a[:,4]cos(θi)^2)

    e2 = a[:,3]cos(θi)^2 + a[:,4]sin(θi)^2
    tj = (sin(2θi) * (a[:,1]-a[:,2])/2)
    tI = sin(2θi)

    #c = ((1+cos(2θi))e1 .+ (1-cos(2θi))e2) ./ tj ./tI

    c = ((1 - cos(2θi))*e1+(1+cos(2θi))*e2)./tj./tI
    storage[i] = 1/mean(log.(abs.(c)))
end

for (i,θi) in enumerate(θ)
    efj = a[:,1]sin(θi)^2 + a[:,2]cos(θi)^2
    epj1 = a[:,3]cos(θi)^2 + a[:,4]sin(θi)^2
    efj1 = a[:,3]sin(θi)^2 + a[:,4]cos(θi)^2
    epj2 = a[:,5]cos(θi)^2 + a[:,6]sin(θi)^2
    tj = (sin(2θi) * (a[:,1]-a[:,2])/2)
    tj1 = (sin(2θi) * (a[:,3]-a[:,4])/2)
    tI = sin(2θi)
    e3 = ((1-cos(2θi))*tj1.^2)./(efj1 .+(1-cos(2θi))/(1+cos(2θi))*epj2)
    #c = ((1+cos(2θi))e1 .+ (1-cos(2θi))e2) ./ tj ./tI

    c = (efj*(1+cos(2θi)) .+ epj1*(1-cos(2θi))-e3)./tj./tI
    storage[i] = 1/mean(log.(abs.(c)))
end


p = plot(θI/pi,ξ_weak_diag, label = "TMM")
scatter!(p, θ/pi,storage, label = "Curve")
gui()
ylims!(0, 3)
xlabel!("θI = θII (π radian)")
ylabel!("ξ")
title!("Weak Disorder Localization")

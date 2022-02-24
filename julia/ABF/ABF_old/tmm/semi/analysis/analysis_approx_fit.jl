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
a = W*(rand(rng,100000,4).-0.5);

storage = Array{Float64}(undef, 100)
θ = LinRange(0,π/2, 100)
for (i,θi) in enumerate(θ)
    tj = sin(2θi)*abs.(a[:,1] .-a[:,2])
    tI = sin(2θi)
    efj = a[:,1]sin(θi)^2 .+ a[:,2]cos(θi)^2
#Approx1
    #original
    c = (1+cos(2θi))*efj./tj/tI
    #simplified
    c = (tI*sin(2θi)*(a[:,2]-a[:,1]) .+ a[:,1]) ./ (2*sin(θi)^2 * abs.(a[:,2] .- a[:,1]))

    storage[i] = 1/mean(log.(abs.(c)))
end
p = plot(θI/pi,ξ_weak_diag, label = "TMM")
scatter!(p, θ/pi,storage, label = "Curve")

xlabel!("θI = θII (π radian)")
ylabel!("ξ")
title!("Weak Disorder Localization")
ylims!(0,3)

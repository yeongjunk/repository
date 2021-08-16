using MAT
using Gadfly
include("xi_abf2.jl")
dir = "/Users/pcs/data/ABF/rawdata/nu2-tm-semi-detangled"
#read all the data and the parameter
ξ_weak, params_weak = ξ_read(dir*"/config-weak-xi_all.mat")

(E, W, θI, θII) =params_expand(params_weak)

data = Array{Float64}(undef, length(θI), 2      )

for i in 1:length(θI)
        data[i,1] = θI[i]
        data[i,2] = ξ_weak[1,1,i,i]
end
matdata = Dict("data" => data)

fn_out = dir*"/python_format_diag.mat"
matwrite(fn_out, matdata, compress=true)
 # write compressed MATLAB file

using MAT
import JSON
include("xi_abf2_div.jl")

cd("/Users/pcs/Data/ABF/Rough")

files = Dict{Int64,String}()
ξ_div = Dict{Int64,Array{Float64,4}}()
params_div = Dict{Int64,Parameters}()

for i in 1:5
    files[i] = ("config"*"_$(i)"*"-xi_all.mat")
    ξ_div[i], params_div[i] = ξ_read(files[i])
end
for i in 2:5
    ξ_div[1][:,params_div[i].W_first:params_div[i].W_last,:,:] = ξ_div[i][:,params_div[i].W_first:params_div[i].W_last,:,:]
end

config["W_first"] = 0
config["W_last"] = 0

config["xi_all"] = ξ_div[1]

config["xi_all"]  = ξ_all                                           #
fn_all = "configmerge"                             # output filename

matwrite(fn_out, config, compress=true)                             # write compressed MATLAB file

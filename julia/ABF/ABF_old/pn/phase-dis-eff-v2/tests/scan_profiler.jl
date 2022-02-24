using JSON
using JLD
include("../library/abf2_pnscan.jl")

configdir = "/Users/pcs/codes/chain/Ladder/Hopping/ABF/pn/phase-dis-eff-v2/config_sample"
config  = JSON.parsefile(configdir)
p = readconfig(config)
BLAS.set_num_threads(1) #Turn off BLAS multi-threading
V, θ = expand_params(p)

## Test singlescan
p_single = full_to_one(p, 6,1)
df = scan_pn(p_single)
@time df  = scan_pn(p_single)
@profiler df  = scan_pn(p_single)

## Test fullscan
@time for i in 1:p.V_num, j in 1:p.θ_num
    p_single = full_to_one(p, i,j)
    df = scan_pn(p_single) #compile
    fn_out = "V_th_$(i)_$(j).jld"
    dict = Dict("W" => p.W, "V" => V[i], "θ" => θ[j], "data" => df)
    JLD.save(fn_out, dict)
end


@time for i in 1:p.V_num, j in 1:p.θ_num
    p_single = full_to_one(p, 1,1)
    df = scan_pn(p_single) #compile
    fn_out = "V_th_$(i)_$(j).jld"
    dict = Dict("W" => p.W, "V" => V[i], "θ" => θ[j], "data" => df)
    JLD.save(fn_out, dict)
end

# Check result
θ, V = expand_params(p)
scatter!(E, pn_mean[:,1,2])
xlims!(-0.7, 0.7)
plot(pn_hist[4,2])

plot(θ, E_bw[:,1])

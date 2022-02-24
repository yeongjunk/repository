using JSON
using Plots
include("../library/abf2_pnscan.jl")

configdir = "/Users/pcs/codes/chain/Ladder/Hopping/ABF/pn/phase-dis-eff/config_sample"
config  = JSON.parsefile(configdir)
p = readconfig(config)
BLAS.set_num_threads(1) #Turn off BLAS multi-threading

smp_num, pn_mean, pn_var, pn_hist, E_bw = scan_pn(pn_fe_ph,p) #compile
@time smp_num, pn_mean, pn_var, pn_hist, E_bw  = scan_pn(pn_sf_ph,p)
@profiler smp_num, pn_mean, pn_var, pn_hist, E_bw = scan_pn(pn_sf_ph,p)


# Check result
θ, V, E = expand_params(p)
scatter!(E, pn_mean[:,1,2])
xlims!(-0.7, 0.7)
plot(pn_hist[4,2])

plot(θ, E_bw[:,1])

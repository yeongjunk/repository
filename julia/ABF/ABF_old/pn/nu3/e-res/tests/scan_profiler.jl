using JSON
include("../library/abf2_pnscan.jl")
print(pwd())
config  = JSON.parsefile("$(pwd())/nu3/pn_e_hist/config")

p = readconfig(config)   # create the Parameters struct from Dict config


BLAS.set_num_threads(1) #Turn off BLAS multi-threading

smp_num, pn_mean, pn_var, pn_hist  = scan_pn(p)


@profiler smp_num, pn_mean, pn_var, pn_hist= scan_pn_diag(p, logscale = false)
@time smp_num, pn_mean, pn_var, pn_hist= scan_pn(p, logscale = false)

plt = plot(pn_hist[1,4])
plot!(plt, )
xlabel!("PN (E = -1, W = 0.1,  θ = 0.25π)")
ylabel!("# Samples")
savefig(plt, "histogram.png")

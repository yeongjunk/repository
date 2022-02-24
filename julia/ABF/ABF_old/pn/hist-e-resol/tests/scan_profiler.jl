include("../library/abf2_pnscan.jl")

p = Parameters(-1., 0.025, 100, 0.01, 10., 2,false, 0., 0.125, 2, 0., 0.125,2, false, 1234, 1000, 200, 4, 20, 30)
BLAS.set_num_threads(1) #Turn off BLAS multi-threading

smp_num, pn_mean, pn_var, pn_hist  = scan_pn_diag(p, logscale = false)
@profiler smp_num, pn_mean, pn_var, pn_hist= scan_pn_diag(p, logscale = false)
@time smp_num, pn_mean, pn_var, pn_hist= scan_pn(p, logscale = false)

plt = plot(pn_hist[1,4])
plot!(plt, )
xlabel!("PN (E = -1, W = 0.1,  θ = 0.25π)")
ylabel!("# Samples")
savefig(plt, "histogram.png")

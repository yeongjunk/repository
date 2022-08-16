using Plots
include("./lib.jl")
params = ScanParameters()
dct = roag_FRP_scan(params)
p = plot(dct["gamma"], dct["R_mean"], yerror = dct["R_ste"])
savefig(p, "test.pdf")

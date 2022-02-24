# compare data with different number of realization(R) at flatband energy
# R = 5000, R = 10000, R = 20000

using MAT
using Plots
using Dates

include("../library/abf2_iprscanlib.jl")

dir = "/Users/pcs/Data/ABF/ipr/E1"
dir1 = "/R5000"
dir2 = "/R10000"
dir3 = "/R20000"
savedir = "/Rtest"
cd(dir*savedir)

file = "/ipr.mat"
vars1 = matread(dir*dir1*file)
vars2 = matread(dir*dir2*file)
vars3 = matread(dir*dir3*file)

vars1["IPR"] = reshape(vars1["IPR"], 1,51,51)
vars2["IPR"] = reshape(vars2["IPR"], 1,51,51)
vars3["IPR"] = reshape(vars3["IPR"], 1,51,51)
ipr = vcat(vars1["IPR"], vars2["IPR"], vars3["IPR"])

W_ipr = collect(LinRange(vars1["W_min"], vars1["W_max"], vars1["W_num"]))
θ_ipr = collect(LinRange(vars1["theta_min"][1], vars1["theta_max"][1], vars1["theta_num"][1]))

i = 10
p = plot()
for j = 3
    scatter!(W_ipr[2:end], ipr[j,i,2:end])
end
display(p)

p = config_read(vars)
print(p)

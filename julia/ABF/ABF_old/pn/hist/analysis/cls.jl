using MAT
using Plots
using Dates

include("../library/abf2_iprscanlib.jl")

dir = "/Users/pcs/data/ABF/ipr/rawdata/cls/er025_n100_r20000"


file = "/ipr.mat"
vars = matread(dir*file)

ipr = vars["IPR"]
norm_a = vars["norm_a"]
norm_b = vars["norm_b"]


W_ipr = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
θ_ipr = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))

plt = plot()
for i = 1:51
    plt = scatter(W_ipr[2:end], ipr[i,2:end])
    display(plt)
    sleep(0.2)
end
display(plt)
i = 31
scatter(W_ipr[2:end],norm_a[i,2:end].^2, legend = false)
scatter!(W_ipr[2:end],norm_b[i,2:end].^2, legend = false)
ylims!(0,1.2)

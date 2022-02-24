#Maximum enhancement

using MAT
using Plots
using Dates

include("../library/abf2_pnscan.jl")

dir = "/Users/pcs/Data/ABF/ipr/rawdata"
dir1 = "/E1/er025_N100_R20000"
dir0 = "/E0/er025_N100_R20000"
savedir = "/Users/pcs/Data/ABF/ipr/analysis/enhance/figures"
cd(savedir)

file = "/ipr.mat"
vars0 = matread(dir*dir0*file)
vars1 = matread(dir*dir1*file)

ipr0 = vars0["IPR"]
ipr1 = vars1["IPR"]

W_ipr0 = collect(LinRange(vars0["W_min"], vars0["W_max"], vars0["W_num"]))
θ_ipr0 = collect(LinRange(vars0["theta_min"][1], vars0["theta_max"][1], vars0["theta_num"][1]))


W_ipr1 = collect(LinRange(vars1["W_min"], vars1["W_max"], vars1["W_num"]))
θ_ipr1 = collect(LinRange(vars1["theta_min"][1], vars1["theta_max"][1], vars1["theta_num"][1]))

ipr_E0_max = Array{Float64}(undef, 51,1)
ipr_E0_idx = Array{Float64}(undef, 51,1)

ipr_E1_weak = Array{Float64}(undef, 51,1)
ipr_E1_max = Array{Float64}(undef, 51,1)
ipr_E1_idx = Array{Float64}(undef, 51,1)

for i = 1:51
    ipr_E0_max[i], ipr_E0_idx[i] = findmax(ipr0[i,:])
    ipr_E1_max[i], ipr_E1_idx[i] = findmax(ipr1[i,:])
    ipr_E1_weak[i] = ipr1[i, 2]
    # savefig(p, "fig"*"$(i)")
end
ipr_zero = 4*(1 .+ cospi.(2θ_ipr0).^2).^(-2)

plt = scatter(θ_ipr0, ipr_E0_max, label = "absolute max", legends =:topleft)
scatter!(plt, θ_ipr0, ipr_E1_weak, label = "weak disorder")
scatter!(plt, θ_ipr0, ipr_zero, label = "CLS")
title!("Absolute MAX of PN in comparison to weak disorder PN")
xlabel!("θ (π rad)")
ylabel!("PN")
savefig(plt, "enhancement")


ipr_rel_max1 = ipr_E0_max./ipr_E1_weak
ipr_rel_max2 = ipr_E0_max./ipr_zero
max_enhance1, enhance_idx1 = findmax(ipr_rel_max1)
max_enhance2, enhance_idx2 = findmax(ipr_rel_max2)

plt = plot(θ_ipr0, ipr_rel_max1)
scatter!(plt, [θ_ipr0[enhance_idx1]], [ipr_rel_max1[enhance_idx1]])
xlabel!("θ (π rad)")
title!("Relative Maximum of PN")
ylabel!("Enhancement")
savefig(plt, "enhancement_ratio")

plt = plot(θ_ipr0, ipr_rel_max2)
scatter!(plt, [θ_ipr0[enhance_idx2]], [ipr_rel_max2[enhance_idx2]])
xlabel!("θ (π rad)")
title!("Relative Maximum of PN (w.r.t. CLS)")
ylabel!("Enhancement w.r.t. CLS")
savefig(plt, "enhancement_ratio_cls")
display(plt)


anim = Plots.Animation()

for i = 1:51
    p = plot()
    scatter!(p, W_ipr0, ipr0[i,:], label = "E = 0")
    scatter!(p, W_ipr1[2:end], ipr1[i,2:end], label = "E = -1")
    scatter!(p, [0], [ipr_zero[i]], label = "CLS")
    display(p)
    sleep(0.2)
    title!("θ=$(θ_ipr0[i])π")
    xlabel!("W")
    ylabel!("PN")
    Plots.frame(anim)
    savefig(p, "fig"*"$(i)")

end

gif(anim, "fig_anim.gif", fps = 10)

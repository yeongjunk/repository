#Weak disorder, effect of changing R and E_window

using MAT
using Plots
using Dates

include("../library/abf2_pnscanlib.jl")

basedir = "/Users/pcs/Data/ABF/ipr/rawdata/e1/"
foldernames = ["er0125_n100_R20000";
"er025_n100_R20000";
"er05_n100_R20000";
]


savedir = "/Users/pcs/Data/ABF/ipr/analysis/var_window"
cd(savedir*"/figures")



labels = ["EW=0.05"; "EW=0.1"; "EW=0.2"]

file = "/ipr.mat"

ipr = Array{Float64}(undef, 51, 51, length(foldernames))
global W = Array{Float64}(undef, 51)
global θ = Array{Float64}(undef, 51)

for i in 1:length(foldernames)
    vars = matread(basedir*foldernames[i]*file)
    ipr[:,:,i] = vars["IPR"]
    if i == 1
        global W = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
        global θ = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))
    end
end

anim = Plots.Animation()
# Plots.frame(anim)
# gif(anim, "fig_anim.gif", fps = 10)

for i = 1:51;
    plt = plot()
    for k = 1:length(foldernames)
        plot!(plt, W[2:end], ipr[i,2:end,k], legends =:bottomright, label=labels[k])
    end
    ipr_zero = 4*(1 .+ cospi.(2θ).^2).^(-2)
    scatter!(plt, [0], [ipr_zero[i]], label = "CLS")
    title!("θ=$(θ[i])π")
    xlabel!("W");ylabel!("PN")
    display(plt)
    Plots.frame(anim)

    savefig(plt, "fig"*"$(i)")
    sleep(0.2)
end
gif(anim, "fig_anim.gif", fps = 10)

#Weak disorder, effect of changing R and E_window

using MAT
using Plots
using Dates

include("../library/abf2_pnscan.jl")

basedir = "/Users/pcs/Data/ABF/ipr/rawdata/e1/"
foldernames = ["er0025_n100_R20000_weak";
"er005_n100_R20000_weak";
"er01_n100_R20000_weak";
"er025_n100_R20000_weak";
"er05_n100_R20000_weak";
]

dir_linscale = "/Users/pcs/Data/ABF/ipr/rawdata/e1/er025_N100_r20000"


savedir = "/Users/pcs/Data/ABF/ipr/analysis/var_window_weak"
cd(savedir*"/figures")



labels = ["EW=0.01"; "EW=0.02"; "EW=0.04"; "EW=0.1"; "EW=0.2"]

file = "/ipr.mat"

ipr = Array{Float64}(undef, 51, 28, length(foldernames))
global W = Array{Float64}(undef, 28)
global θ = Array{Float64}(undef, 51)

for i in 1:length(foldernames)
    vars = matread(basedir*foldernames[i]*file)
    ipr[:,:,i] = vars["IPR"]
    if i == 1
        global W = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
        global θ = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))
    end
end

vars = matread(dir_linscale*file)
ipr_lin = vars["IPR"]
W_lin = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
θ_lin = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))
ipr_zero = 4*(1 .+ cospi.(2θ).^2).^(-2)

anim = Plots.Animation()
for i = 1:51;
    plt = plot()
    for k = 1:length(foldernames)
        plot!(plt, 10 .^W[1:end], ipr[i,1:end,k], legends =:bottomright, label=labels[k])
    end
    plot!(plt, W_lin[2:end], ipr_lin[i,2:end] , label = "full scan(EW=0.1)", color=:black)

    scatter!(plt, [0], [ipr_zero[i]], label = "CLS")
    title!("θ=$(θ[i])")
    xlabel!("W");ylabel!("PN")
    xlims!(-0.1,4)
    display(plt)
    sleep(0.2)
    Plots.frame(anim)
end
gif(anim, "weak_W_linscale.gif", fps = 10)


anim = Plots.Animation()

for i = 1:51
    plt = plot()
    # scatter!(plt, [-5], [ipr_zero[i]], label = "CLS")
    for k = 1:length(foldernames)
        plot!(plt, W[1:end], ipr[i,1:end,k], legends =:topleft, label=labels[k])
        title!("θ=$(θ[i])")
        xlabel!("log10W");ylabel!("PN")
    end
    display(plt)
    sleep(0.1)
    Plots.frame(anim)
    savefig(plt, "fig_logscale"*"$(i)")
end
gif(anim, "weak_W_logscale.gif", fps = 10)

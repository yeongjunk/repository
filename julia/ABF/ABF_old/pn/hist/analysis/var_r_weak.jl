#Weak disorder, effect of changing R and E_window

using MAT
using Plots
using Dates

include("../library/abf2_iprscanlib.jl")

basedir = "/Users/pcs/Data/ABF/ipr/rawdata/e1/"
foldernames = ["er025_n100_R100_weak" "er025_n100_R1000_weak"
 "er025_n100_R10000_weak" "er025_n100_R20000_weak"]

realizations = ["R=100" "R=1000" "R=10000" "R=20000"]

file = "/ipr.mat"

ipr = Array{Float64}(undef, 51, 28, length(foldernames))
for i in 1:length(foldernames)
    vars = matread(basedir*foldernames[i]*file)
    ipr[:,:,i] = vars["IPR"]
    if i == 1
        W = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
        θ = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))
    end
end

savedir = "/Users/pcs/Data/ABF/ipr/analysis/weak"
mkdir(savedir*"/figures")
cd(savedir*"/figures")


i = 51;
plt = plot()
for k = 1:length(foldernames)
    scatter!(plt, 10 .^W[1:end], ipr[i,1:end,k], legends =:topleft, label=realizations[k])
end
display(plt)
for i = 1:51
    plt = plot()
    scatter!(plt, W_ipr, ipr[i,:], legends = false)
    title!("θ=$(θ_ipr[i])π")
    xlabel!("W")
    ylabel!("PN")
    savefig(plt, "fig"*"$(i)")
end

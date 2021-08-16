#Maximum enhancement

using JLD
using Plots
using Dates

include("../library/abf2_pnscan.jl")
dir = [
"/Users/pcs/Data/ABF/pn/rawdata/hist/e1_er025_weak/",
"/Users/pcs/Data/ABF/pn/rawdata/hist/e1_er0025_weak/",
"/Users/pcs/Data/ABF/pn/rawdata/hist/e1_er00025_weak/"]

file = "pn.jld"
savedir = "/Users/pcs/Data/ABF/pn/analysis/histogram/weak"
cd(savedir)

vars = [load(dir[1]*file),
load(dir[2]*file),
load(dir[3]*file)]

W = collect(LinRange(vars[1]["W_min"], vars[1]["W_max"], vars[1]["W_num"]))
θ = collect(LinRange(vars[1]["theta_min"][1], vars[1]["theta_max"][1], vars[1]["theta_num"][1]))

pn_mean = [vars[1]["pn_mean"], vars[2]["pn_mean"], vars[3]["pn_mean"]]
smp_num = [vars[1]["smp_num"], vars[2]["smp_num"], vars[3]["smp_num"]]
pn_var = [vars[1]["pn_var"], vars[2]["pn_var"], vars[3]["pn_var"]]
pn_hist = [vars[1]["pn_hist"], vars[2]["pn_hist"], vars[3]["pn_hist"]]


pn_cls = 4*(1 .+ cospi.(2θ).^2).^(-2)
#
# if isdir("std") == false
#     mkdir("stdn")
# end
anim = Plots.Animation()
for i = 1:length(θ)
    p = plot()
    scatter!(p, W, sqrt.(pn_var[1][i,:]), label = "E = -1")
    scatter!(p, W, sqrt.(pn_var[2][i,:]), label = "E = -1")
    scatter!(p, W, sqrt.(pn_var[3][i,:]), label = "E = -1")
    title!("θ=$(θ[i])π")
    xlabel!("W")
    ylabel!("σ_PN")
    display(p)
    sleep(0.2)
    Plots.frame(anim)
    # savefig(p, "std/fig_std"*"$(i)")
end
# gif(anim, "std/fig_anim.gif", fps = 10)

anim = Plots.Animation()
for i = 1:length(θ)
    p = plot()
    scatter!(p, W, (pn_mean[1][i,:]), label = "E = -1")
    scatter!(p, W, (pn_mean[2][i,:]), label = "E = -1")
    scatter!(p, W, (pn_mean[3][i,:]), label = "E = -1")
    title!("θ=$(θ[i])π")
    xlabel!("W")
    ylabel!("σ_PN")
    display(p)
    sleep(0.2)
    Plots.frame(anim)
    # savefig(p, "std/fig_std"*"$(i)")
end

anim = Plots.Animation()
for i = 1:length(θ)
    p = plot()
    scatter!(p, W, pn_mean[1][i,:], label = "E = -1", yerr = sqrt.(pn_var[1][i,:]))
    scatter!(p, W, pn_mean[2][i,:], label = "E = -1", yerr = sqrt.(pn_var[2][i,:]))
    scatter!(p, W, pn_mean[3][i,:], label = "E = -1", yerr =sqrt.(pn_var[3][i,:]))
    title!("θ=$(θ[i])π")
    xlabel!("W")
    ylabel!("σ_PN")
    display(p)
    sleep(0.2)
    Plots.frame(anim)
    # savefig(p, "std/fig_std"*"$(i)")
end


#
if isdir("gifs") == false
    mkdir("gifs")
end
for j = 1:2
    anim1 = Plots.Animation()
    anim2 = Plots.Animation()
    anim3 = Plots.Animation()

    for k = 1:28
        p1 = plot(pn_hist[1][j,k], color =:black, label = "ER025")
        vline!([pn_cls[j]], color =:red, label = "CLS")
        vline!([pn_mean[1][j,k]], color =:blue, label = "Mean" )
        xlabel!("PN")
        ylabel!("Number of samples")
        title!("θ = $(θ[j]), W = $(round(10^W[k],digits = 4))")
        display(p1)
        sleep(0.2)
        Plots.frame(anim1)

        p2 = plot(pn_hist[2][j,k], color =:black, label = "ER0025")
        vline!([pn_cls[j]], color =:red, label = "CLS")
        vline!([pn_mean[2][j,k]], color =:blue, label = "Mean" )
        xlabel!("PN")
        ylabel!("Number of samples")
        title!("θ = $(θ[j]), W = $(round(10^W[k],digits = 4))")
        display(p2)
        sleep(0.2)
        Plots.frame(anim2)

        p3 = plot(pn_hist[3][j,k], color =:black, label = "ER00025")
        vline!([pn_cls[j]], color =:red, label = "CLS")
        vline!([pn_mean[3][j,k]], color =:blue, label = "Mean" )
        xlabel!("PN")
        ylabel!("Number of samples")
        title!("θ = $(θ[j]), W = $(round(10^W[k],digits = 4))")
        display(p3)
        sleep(0.01)
        Plots.frame(anim3)
    end
    gif(anim1, "gifs/fig_anim1_$(j).gif", fps = 10)
    gif(anim2, "gifs/fig_anim2_$(j).gif", fps = 10)
    gif(anim3, "gifs/fig_anim3_$(j).gif", fps = 10)
end



y_ran = midpoints(pn_hist[51,2].edges[1])

data = Array{Float64}(undef, length(y_ran),50)
j = 20
for i in 1:50
    data[:,i] = pn_hist[j,i+1].weights/norm(pn_hist[j,i+1].weights)
end

heatmap(W[2:end],y_ran[1:250], data[1:250,:])

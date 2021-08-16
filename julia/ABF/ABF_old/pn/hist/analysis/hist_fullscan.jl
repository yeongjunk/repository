#Maximum enhancement

using JLD
using Plots
using Dates

include("../../library/abf2_pnscan.jl")
dir = "/Users/pcs/Data/ABF/pn/rawdata/hist/e1_er025/1/"
file = "pn.jld"
savedir = "/Users/pcs/Data/ABF/pn/analysis/histogram/"
cd(savedir)

vars = load(dir*file)

println(vars["info"])
W = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
θ = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))

pn_mean = vars["pn_mean"]
smp_num = vars["smp_num"]
pn_var = vars["pn_var"]
pn_hist = vars["pn_hist"]

pn_cls = 4*(1 .+ cospi.(2θ).^2).^(-2)

if isdir("std") == false
    mkdir("stdn")
end
anim = Plots.Animation()
for i = 1:length(θ)
    p = plot()
    scatter!(p, W[2:end], sqrt.(pn_var[i,2:end]), label = "E = -1")
    display(p)
    sleep(0.2)
    title!("θ=$(θ[i])π")
    xlabel!("W")
    ylabel!("σ_PN")
    Plots.frame(anim)
    savefig(p, "std/fig_std"*"$(i)")
end
gif(anim, "std/fig_anim.gif", fps = 10)

if isdir("mean") == false
    mkdir("mean")
end
anim = Plots.Animation()
for i = 1:length(θ)
    p = plot()
    scatter!(p, W[2:end], pn_mean[i,2:end], label = "E = -1", color = :red)
    scatter!(p,[0], [pn_cls[i]], label = "CLS", color = :green)
    display(p)
    sleep(0.1)
    title!("θ=$(θ[i])π")
    xlabel!("W")
    ylabel!("mean(PN)")
    Plots.frame(anim)
    savefig(p, "mean/fig_mean"*"$(i)")
end
gif(anim, "mean/fig_anim.gif", fps = 10)

if isdir("mean_errbar") == false
    mkdir("mean_errbar")
end
anim = Plots.Animation()
for i = 1:length(θ)
    p = plot()
    plot!(p, W[2:end], (pn_mean[i,2:end]), yerr = sqrt.(pn_var[i,2:end]), label = "E = -1", color =:red)
    scatter!(p,[0], [pn_cls[i]], label = "CLS", color = :green)
    sleep(0.2)
    title!("θ=$(θ[i])π")
    xlabel!("W")
    ylabel!("mean(PN)")
    display(p)
    Plots.frame(anim)
    savefig(p, "mean_errbar/fig_mean_errbar"*"$(i)")
end
gif(anim, "mean_errbar/fig_anim.gif", fps = 10)

if isdir(gifs) == false
    mkdir(gifs)
end
for j = 1:51
    anim = Plots.Animation()
    for k = 1:51
        ymax,_ = findmax(pn_hist[j,2].weights)
        p = plot(pn_hist[j,k], color =:black, label = "PN histogram")
        vline!([pn_cls[j]], color =:red, label = "CLS")
        vline!([pn_mean[j,k]], color =:blue, label = "Mean" )
        xlabel!("PN")
        ylabel!("Number of samples")
        ylims!(0, ymax)
        title!("θ = $(θ[j]), W = $(round(W[k],digits = 2))")
        # display(p)
        Plots.frame(anim)
    end
    gif(anim, "gifs/fig_anim$(j).gif", fps = 10)
end



y_ran = midpoints(pn_hist[51,2].edges[1])

data = Array{Float64}(undef, length(y_ran),50)
j = 20
for i in 1:50
    data[:,i] = pn_hist[j,i+1].weights/norm(pn_hist[j,i+1].weights)
end

heatmap(W[2:end],y_ran[1:250], data[1:250,:])

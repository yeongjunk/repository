#Maximum enhancement

using JLD
using Plots
using Dates
using Polynomials
include("../library/abf2_pnscan.jl")

file = "/Users/pcs/Data/ABF/pn/rawdata/hist_e/fullscan/1/pn.jld"
vars = load(file)
PN_mean = vars["PN_mean"]
N = vars["N"]
PN_var = vars["PN_var"]
PN_hist = vars["PN_hist"]
θ = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))
for i in 2:10
    file = "/Users/pcs/Data/ABF/pn/rawdata/hist_e/fullscan/$(i)/pn.jld"
    vars = load(file)
    tPN_mean = vars["PN_mean"]
    tN = vars["N"]
    tPN_var = vars["PN_var"]
    tPN_hist = vars["PN_hist"]
    tθ = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))


    PN_mean = vcat(PN_mean, tPN_mean)
    N = vcat(N, tN)
    PN_var = vcat(PN_var, tPN_var)
    PN_hist = vcat(PN_hist, tPN_hist)
    θ = vcat(θ, tθ)
end

PN_cls = 4*(1 .+ cospi.(2θ).^2).^(-2)



println(vars["info"])
W = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
E_bin_num = vars["E_bin_num"]
E_bins = range(-0.5, 0.5, length = E_bin_num)


PN_mean_E0 = Array{Float64}(undef, 50, 28)
for j in 1:length(θ)
    for i in 1:length(W)
        E = (E_bins * W[i]) .-1
        x = findall(x -> x > 0., E)
        if x[end] == length(E)
            PN_mean_E0[j,i] = NaN
        else
            PN_mean_E0[j,i] = PN_mean[j,i,x[1]]
        end
    end
end

for i = 1:length(θ)
    p = plot()
    plot!(p, W, (PN_mean[i,:,15]), yerr = sqrt.(PN_var[i,:,15]),
    label = "E = -1", color =:red,
    linewidth = 3,
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize=18,
    yguidefontsize=18,
    titlefontsize= 18,
    legendfontsize=10,
    framestyle = :box
    )
    scatter!([0], [PN_cls[i]], color =:blue, label = "CLS")
    ylabel!("<PN>")
    xlabel!("W")
    display(p)
    #savefig(p, "/Users/pcs/Data/ABF/analysis/pn/"*"$(i)")
end


for i = 1:length(θ)
    p = plot()

    plot!(p, W, PN_mean_E0[i,:],
    label = "E = 0", color =:red,
    linewidth = 3,
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize=18,
    yguidefontsize=18,
    titlefontsize= 18,
    legendfontsize=10,
    framestyle = :box
    )
    plot!(p, W, PN_mean[i,:,15], label = "E = -1",
    color = :blue,
    linewidth = 3,
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize=18,
    yguidefontsize=18,
    titlefontsize= 18,
    legendfontsize=10,
    framestyle = :box

    )
    scatter!([0], [PN_cls[i]], color =:green, label = "CLS")
    ylabel!("<PN>")
    xlabel!("W")
    display(p)
    savefig(p, "/Users/pcs/Data/ABF/analysis/pn/"*"$(i)")
end

PN_max = Array{Float64}(undef, length(θ))
PN_mean_E0 = map(x -> isnan(x) ? zero(x) : x, PN_mean_E0)
for i = 1:length(θ)
    (PN_max[i], _) = findmax(PN_mean_E0[i,:])
end

p = plot(θ, PN_max, label = "PN_max",
color =:red,
linewidth = 3,
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
legendfontsize=10,
framestyle = :box,
legends = :topleft
)
xlabel!("θ")
ylabel!("PN")
plot!(p,θ, PN_mean[:,1,15], label = "PN_weak",
color =:blue,
linewidth = 3)

plot!(p,θ, PN_cls, label = "PN_cls",
color = :green,
linewidth = 3)

savefig(p, "/Users/pcs/Data/ABF/analysis/pn/"*"comparison")

p = plot(θ, PN_max./PN_mean[:,1,15], label = "Enhancement",
color =:red,
linewidth = 3,
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
legendfontsize=10,
framestyle = :box,
legend = :topright
)
xlabel!("θ")
ylabel!("Enhancement")
maxratio, maxidx = findmax(PN_max./PN_mean[:,1,15])
println(θ[maxidx])

savefig(p, "/Users/pcs/Data/ABF/analysis/pn/"*"enhancement")



E = (E_bins * W[1]) .-1
p = plot(E, PN_mean[50,1,:],
color =:red,
linewidth = 3,
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
legendfontsize=10,
framestyle = :box,
legend = :topright)
ylabel!("<PN>")
xlabel!("E")

savefig(p, "/Users/pcs/Data/ABF/analysis/pn/"*"PNE")

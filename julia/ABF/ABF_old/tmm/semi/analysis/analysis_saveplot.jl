using MAT
using Random
using Plots
using Dates
using Distributions
using LaTeXStrings
include("xi_abf2.jl")

# Go to the directory where the data is
file_full = "/Users/pcs/Data/ABF/tmm/fullscan/config-xi_all.mat"
file_weak = "/Users/pcs/Data/ABF/tmm/fullscan/config-weak-xi_all.mat"

#read all the data and the parameter
ξ, params = ξ_read(file_full)
ξ_weak, params_weak = ξ_read(file_weak)
# create paraameter space from the configuration
(E, W, θI, θII) =params_expand(params)

ξ_max = Array{Float64}(undef, params.θI_num,params.θII_num)
ξ_maxidx = Array{CartesianIndex}(undef, params.θI_num,params.θII_num)
for i in 1:params.E_num, j in 1:params.E_num
    ξ_max[i,j] ,ξ_maxidx[i,j] = findmax(ξ[:,:,i,j])
end

ξ_weak_diag = [ξ_weak[1,1,i,i] for i in 1:params.θI_num]
ξ_max_diag = [ξ_max[i,i] for i in 1:params.θI_num]
ξ_relmax = ξ_max./ξ_weak[1,1,:,:]


# weak disorder, ξ(θI,θII)
p_weak = heatmap();
heatmap!(p_weak,θII/π  ,θI/π, ξ_weak[1,1,:,:],
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
)
title!("ξ, at W = 0.01, E = -1");xlabel!("θI (π rad)");ylabel!("θII (π rad)")
savefig(p_weak, "/Users/pcs/Data/ABF/analysis/tmm/weak")

# weak disorder, ξ(θ) (diagonal)
p_weak_diag = plot();
plot!(p_weak_diag, θII/π, ξ_weak_diag,
linewidth = 3, legend = false,
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
legendfontsize=10,
framestyle = :box
)
xlabel!("θ (π rad)");ylabel!("ξ")
savefig(p_weak_diag, "/Users/pcs/Data/ABF/analysis/tmm/weak_diag")





θ_idx = [2 8 14 20 26];
#ξ(E) plot.
W_idx = [1 6 11 16 21 26 31];
for i in 1:length(θ_idx)
    for k in 1:length(θ_idx)
        p_E = plot()
        for j = 1:length(W_idx)
            plot!(p_E,E, ξ[:,W_idx[j],θ_idx[i],θ_idx[k]],
            label="W = $(round(W[W_idx[j]],digits = 2))",
            legend =:bottomright,
            linewidth = 3,
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize=18,
            yguidefontsize=18,
            titlefontsize= 18,
            legendfontsize=10,
            framestyle = :box
            )
        end
    angles  = "(θI, θII) = ($(round(θI[θ_idx[i]]/π,digits = 2))π, $(round(θII[θ_idx[k]]/π,digits = 2))π)"
    annotate!((-0.6, Plots.ylims(p_E)[1] + 0.05*(Plots.ylims(p_E)[2]-Plots.ylims(p_E)[1]), angles))
    xlabel!("E, Energy")
    ylabel!("ξ, Localization Length")
    display(p_E)
    savefig(p_E, "/Users/pcs/Data/ABF/analysis/tmm/xi_E/fig$(i)_$(k)")
    end
end



#ξ(E) plot.(Fig 2) index: i for θ, j for W. 5 plots(depending on W) in a same figure.
E_idx = [1 6 11 16 21 26];
for i in 1:length(θ_idx), k in 1:length(θ_idx)
    p_W = plot()
    for j = 1:length(E_idx)
        plot!(p_W,W, ξ[E_idx[j],:,θ_idx[i],θ_idx[k]],
        label="E = $(round(E[E_idx[j]],digits = 2))",
        legend =:topright,
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize=18,
        yguidefontsize=18,
        titlefontsize= 18,
        legendfontsize=10,
        framestyle = :box
        )
    end
    angles  = "(θI, θII) = ($(round(θI[θ_idx[i]]/π,digits = 2))π, $(round(θII[θ_idx[k]]/π,digits = 2))π)"
    xlabel!("W, Disorder strength")
    ylabel!("ξ, Localization length")
    annotate!((6, Plots.ylims(p_W)[1] + 0.05*(Plots.ylims(p_W)[2]-Plots.ylims(p_W)[1]), angles))
    display(p_W)
    savefig(p_W, "/Users/pcs/Data/ABF/analysis/tmm/xi_W/fig$(i)_$(k)")
end


#ξ(E,W) full scan
for i in 1:length(θI)
    anim = Plots.Animation()
    for j in 1:length(θII)
        p_EW = heatmap(W,E,ξ[:,:,i,j],
        title = "θI= $(round(θI[i]/π,digits = 2))π, θII= $(round(θII[j]/π,digits = 2))π",
        xlabel = "W",
        ylabel = "E",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize=18,
        yguidefontsize=18,
        titlefontsize= 18,
        )
        display(p_EW)
        Plots.frame(anim)
    end
    gif(anim, "/Users/pcs/Data/ABF/analysis/tmm/xi_EW/fig$(i).gif", fps = 10)
end


#ξ(E,W) diagonal angle
anim = Plots.Animation()
for i in 1:length(θI)
    p_EW = heatmap(W,E,ξ[:,:,i,i],
    title = "θI= $(round(θI[i]/π,digits = 2))π, θII= $(round(θII[i]/π,digits = 2))π",
    xlabel = "W",
    ylabel = "E",
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize=18,
    yguidefontsize=18,
    titlefontsize= 18,
    )
    display(p_EW)
    Plots.frame(anim)
end
gif(anim, "/Users/pcs/Data/ABF/analysis/tmm/fullscan_result/xi_EW/fig_diag.gif", fps = 10)


# Absolute maximum
p_EW_abs_max = heatmap(θII/pi,θI/pi, ξ_max,
linewidth = 3, legend = false,
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
legendfontsize=10,
framestyle = :box
)
xlabel!("θI (π radian)");ylabel!("θII (π radian)");title!("Absolute Maximum of ξ")
savefig(p_EW_abs_max, "/Users/pcs/Data/ABF/analysis/tmm/absmax_2D.png")

# Absolute max, (diagonal)
p_max_diag = plot();
plot!(p_max_diag, θII/π, ξ_max_diag,
linewidth = 3, legend = false,
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
legendfontsize=10,
framestyle = :box
)
xlabel!("θ (π rad)");ylabel!("ξ")
savefig(p_max_diag, "/Users/pcs/Data/ABF/analysis/tmm/max_diag")



# Enhancement
p_EW_rel_max = heatmap(θII/pi,θI/pi, ξ_relmax,
linewidth = 3, legend = false,
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
legendfontsize=10,
framestyle = :box
)

xlabel!("θI (π radian)");ylabel!("θII (π radian)");title!("Enhancement")
savefig(p_EW_rel_max, "/Users/pcs/Data/ABF/analysis/tmm/enhancement_2D.png")

relmaxval, relmaxidx = findmax(ξ_relmax[1:div(end,2),1:div(end,2)])

println("relative maximum at θI = θII = $(θII[relmaxidx[2]]/π)π")
println("relative maximum value $(relmaxval)")

ξ_diagrelmax = Array{Float64}(undef, params.θI_num)

for i in 1:params.θI_num
    ξ_diagrelmax[i] = ξ_relmax[i,i]
end
p_diag_rel_max = plot(θI/pi, ξ_diagrelmax,
linewidth = 3, legend = false,
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
legendfontsize=10,
framestyle = :box
)
xlabel!("θ (π radian)")
ylabel!("ξ")
savefig(p_diag_rel_max, "/Users/pcs/Data/ABF/analysis/tmm/enhancement.png")


p_Wmaxidx = heatmap(θII/π, θI/π, W[getindex.(ξ_maxidx[:,:],2)])
xlabel!("θI (π radian)");ylabel!("θII (π radian)");title!("Value of W where ξ is maximum")
println("Minimum: $(findmin(W[getindex.(ξ_maxidx[:,:],2)]))")
println("Maximum: $(findmax(W[getindex.(ξ_maxidx[:,:],2)]))")




pt = plot()
plot!(pt,θI/pi,ξ_weak_diag, label = "weak disorder", linewidth = 3,
)
plot!(pt,θI/pi,ξ_max_diag, label = "maximum",
linewidth = 3,
xtickfontsize = 14,
ytickfontsize = 14,
xguidefontsize=18,
yguidefontsize=18,
titlefontsize= 18,
legendfontsize=10,
framestyle = :box
)
xlabel!("θ=θI=θII (π rad)")
ylabel!("ξ")
savefig(pt, "/Users/pcs/Data/ABF/analysis/tmm/comparison.png")


pt3 = heatmap()
heatmap!(pt3, θI/pi, θI/pi, ξ_max - ξ_weak[1,1,:,:])

cd(foldername)
#Save figures
#
# savefig(p_W,"fig2")
# savefig(pEW,"fig3")
# savefig(p_EW_abs_max, "fig4")
# savefig(p_EW_rel_max, "fig5")
# savefig(p_Wmaxidx, "fig6")
# savefig(p_diag_rel_max, "fig7")
# savefig(pt, "fig8")
# savefig(pt2, "fig9")
# savefig(pt3, "fig10")

using Interpolations
using LaTeXStrings
using Plots
using MAT
using StatsBase
using Statistics
using Plots.PlotMeasures

function relative(f, r; sp)
    p = plot!()
    lims = f(p[sp])
    return lims[1] + r * (lims[2]-lims[1])
end
relativex(r; sp::Int=1) = relative(Plots.xlims, r; sp=sp)
relativey(r; sp::Int=1) = relative(Plots.ylims, r; sp=sp)


gr(dpi = 300)
savedir = "/Users/pcs/data/ABF3D/analyzed-data/"

# ROAG parameters
L = [10 15 20 25 30]
R = [2000 800 250 100 60]
θ = range(0.0001, 0.25, length = 26)
E_edges = collect(0.65:0.1:1.35)
pushfirst!(E_edges,0.625)
pushfirst!(E_edges, 0.6)
push!(E_edges,  1.375)
push!(E_edges,  1.4)
push!(E_edges, 1.30)
push!(E_edges, 0.70)
sort!(E_edges)
## Roag_full Parameters
L_m = [10 15 20 25 30]
R_m = [2000 800 250 100 60]
θ_m = range(0.0001, 0.25, length = 26)

E_edges_m = collect(0.65:0.1:1.35)
pushfirst!(E_edges_m,0.625)
pushfirst!(E_edges_m, 0.6)
push!(E_edges_m,  1.375)
push!(E_edges_m,  1.4)
push!(E_edges_m, 1.30)
push!(E_edges_m, 0.70)
sort!(E_edges_m)
## ROAG_fine Parameters
L_f = [10 15 20 25]
R_f = [10000 4000 600 200]
θ_f = range(0.09, 0.11, length = 11)
E_edges_f = collect(range(0.95, 1.05, length = 2))

## PN_full parameters
L_pn = [10 15 20 25 30]
R_pn = [1200 400 100 65 50]
dir_pn = "/Users/pcs/data/ABF3D/full-ipr/"
θ_pn = range(0.0001, 0.25, length = 26)
savedir_pn = "/Users/pcs/data/ABF3D/analyzed-data/"
E_edges_pn = collect(range(0.6, 1.4, length = 20))
# E_edges_pn = E_edges_m
## Data
r_m_dir = "/Users/pcs/data/ABF3D/analyzed-data/roag_merged.mat"
pn_dir = "/Users/pcs/data/ABF3D/analyzed-data/pn_old.mat"

r_m_data = matread(r_m_dir)
pn_data = matread(pn_dir)
r_m_m = r_m_data["r_m_m"]
r_m_m_f = r_m_data["r_m_m_f"]

θ_merged = r_m_data["th_merged"]
r_m_m_merged = r_m_data["r_m_m_merged"]
pn_exp_diff = pn_data["pn_exp_diff"]
## Plot settings
marker = (:circle, 4, 1., stroke(-0.5, 1., :black))
line = (:line, :solid, 2)
palette_roag = :Dark2_5
default(
    framestyle = :box,
    size = (600,400),
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    legend = (0.1, 0.75),
    fontfamily = "computer modern",
    tickfontsize = 18,
    guidefontsize = 22,
    legendfontsize = 16,
    )
## Figure 3 (a)
p_a = plot(
    xlabel = L"$\theta / \pi$",
    ylabel = L"$E$",
    xlims = (-0.005, 0.255),
    ylims = (0.6, 1.4),
    xticks = range(0.,0.25, step = 0.05),
    yticks = range(0.6, 1.4, length = 5)
)
# Interpolations
itp = interpolate(pn_exp_diff[:, :, end-2], BSpline(Linear()))
r_new = itp(range(1, length(E_edges_pn)-1, length = 100), range(1,26, length = 100))
θ_itp = interpolate(θ_pn, BSpline(Linear()))
θ_new  = θ_itp(range(1,length(θ_pn), length = 100))
E_itp = interpolate(midpoints(E_edges_pn), BSpline(Linear()))
E_new = E_itp(range(1,length(E_itp), length = 100))


heatmap!(p_a, θ_new, E_new, r_new, c= :delta)
savefig(p_a, savedir*"alpha_full.png")
plot!(colorbar = false)
savefig(p_a, savedir*"alpha_full_nocbar.png")
p_a


## Figure 3 (b)
itp = interpolate(r_m_m[:,:,end], BSpline(Linear()))
r_new = itp(range(1, length(E_edges_m)-1, length = 100), range(1,26, length = 100))
θ_itp = interpolate(θ, BSpline(Linear()))
θ_new  = θ_itp(range(1,length(θ), length = 100))
E_itp = interpolate(midpoints(E_edges_m), BSpline(Linear()))
E_new = E_itp(range(1,length(E_itp), length = 100))

p_b = plot(
    xlabel = L"$\theta / \pi$",
    ylabel = L"$E$",
    xlims = (-0.005, 0.255),
    ylims = (0.6, 1.4),
    xticks = range(0.,0.25, step = 0.05),
    yticks = range(0.6, 1.4, length = 5)
)
heatmap!(p_b, θ_new, E_new, r_new, c= :balance)
savefig(p_b, savedir*"roag_full.png")
plot!(colorbar = false)
savefig(p_b, savedir*"roag_full_nocbar.png")


## Figure 3 (c)
k = 10
E_bins_k = round.(midpoints(E_edges_pn)[k], digits = 2)

p_c = plot(
    xlabel =  L"$\theta / \pi$",
    xlims = (0., 0.25),
    xticks = 0.0:0.05:0.25,
    ylabel =  L"$\alpha$",
    ylims = (-0.2, 3),
    yticks = 0.0:1.0:3.0
)
# annotate!(0.055, 1.2, Plots.text(L"E = %$(E_bins_k)",18, :right, :black))

for i in 1:size(pn_exp_diff, 3)
    plot!(p_c, θ_pn, pn_exp_diff[k,:,i],
        marker = marker,
        line = line, c = i,
        label = L"L: %$(L[i]) \rightarrow %$(L[i+1])")
end
for i in 1:size(pn_exp_diff, 3)
    scatter!(p_c, θ_pn, pn_exp_diff[k,:,i],
        marker = marker,
        c = i,
        label = :none)
end
# vline!(p,[0.1], c = :black, lw = 2, linestyle = :dot, label = :none)
display(p_c)
savefig(p_c, savedir*"alpha_E$(E_bins_k).pdf")


## Figure 3(d)
k = 8
E_bins_k = round.(midpoints(E_edges)[k], digits = 2)

p_d = plot(palette = palette_roag,
    xlabel =  L"$\theta / \pi$",
    ylabel =  L"$\langle r \rangle$",
    xlims = (0., 0.25),
    ylims = (0.38, 0.54),
    xticks = 0.0:0.05:0.25,
    yticks = range(0.38,0.54, length = 5)
)
# annotate!(0.055, 0.435, Plots.text(L"E = %$(E_bins_k)",18, :right, :black))

for i in 1:4
    plot!(p_d, θ_merged, r_m_m_merged[:,i],
        marker = marker,
        line = line,
        c = i,
        label = L"L = %$(L[i])")
end
for i in 1:4
    scatter!(p_d, θ_merged, r_m_m_merged[:,i],
        marker = marker,
        c = i,
        label = :none)
end

## Figure 3 (d) inset
k = 1
plot!(p_d, sp = 2, palette = palette_roag, legend = false,
    inset = (1, bbox(0.049, 0.25, 0.4, 0.55, :bottom, :right)),
    xlabel =  L"$\theta / \pi$",
    xlims = (0.09, 0.11),
    ylims = (0.49, 0.532),
    xticks = 0.09:0.01:0.11,
    yticks = range(0.49, 0.53, length = 3)
)

for i in 1:4
    plot!(p_d, sp = 2, θ_f, r_m_m_f[k,:,i],
        marker = marker,
        line = line,
        c = i,
        label = L"L = %$(L[i])")
end
for i in 1:4
    scatter!(p_d, sp = 2, θ_f, r_m_m_f[k,:,i],
        marker = marker,
        c = i,
        label = :none)
end
    # vline!(p,[0.1], c = :black, lw = 2, linestyle = :dot, label = :none)
vline!(p_d, sp = 2, [0.1], c = :black, lw = 2, linestyle = :dot, label = :none)
display(p_d)
savefig(p_d, savedir*"roag_w_inset.pdf")

println("roag plotted")


p_full = plot(p_a, p_b, p_c, p_d, layout = 4, size = (1200, 800))

annotate!(sp=1,[(relativex(0.03; sp=1), relativey(0.93; sp=1), text("(a)",:white, :left, 22, "computer modern"))])
annotate!(sp=2,[(relativex(0.03; sp=2), relativey(0.93; sp=2), text("(b)",:white ,:left, 22, "computer modern"))])
annotate!(sp=3,[(relativex(0.03; sp=3), relativey(0.93; sp=3), text("(c)",:black ,:left, 22, "computer modern"))])
annotate!(sp=4,[(relativex(0.03; sp=4), relativey(0.93; sp=4), text("(d)",:black ,:left, 22, "computer modern"))])


plot!(p_d,
    left_margin = [4.5mm 0mm],
    right_margin = [4.5mm 0mm],
    )

plot!(p_a,
    left_margin = [5mm 0mm],
    )



savefig(p_full, savedir*"full_3d.pdf")


p_full

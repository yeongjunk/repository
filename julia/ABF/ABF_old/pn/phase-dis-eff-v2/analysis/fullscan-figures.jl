using JSON
using ColorSchemes
using Interpolations
include("../library/abf2_analysis_tool.jl")
include("../library/abf2_pnscan.jl")

# File names and directories
fullscan_config_fn = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-abf-sf/config"
fullscan_dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-abf-sf/R1/"
del0_dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-abf-sf-del0/rawdata/"
delinf_dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-abf-sf-delinf/rawdata/"
savedir = "/Users/pcs/data/ABF/analysis/nu2-pn-phase-dis-eff-new/"

## Load parameters
p_dict = JSON.parsefile(fullscan_config_fn)
p = readconfig(p_dict)
## Load all list of JLD files from dirs
fns_full = readdir_jld(fullscan_dir, join = true)
fns_del0 = readdir_jld(del0_dir, join = true)
fns_delinf = readdir_jld(delinf_dir, join = true)
## Setup E resolutions
binnum = 51
E_edges = range(-0.5, 0.5, length = binnum + 1)
## Load fullscan summary
df_full = summary_full(fns_full, E_edges)
df_del0 = summary_full(fns_del0, E_edges)
df_delinf = summary_full(fns_delinf, E_edges)
df_full.V .= log10.(df_full.V)
## Find maximum of PN(E) from full scan(non-zero finite delta)
V, θ = expand_params(p)
V = log10.(V)
dfmax = DataFrame(E = Float64[], V = Float64[], θ = Float64[], PN_mean = Float64[], PN_err = Float64[], N_E = Int64[])
for i in 1:length(V), j in 1:length(θ)
    df_single = df_full[(df_full.V .== V[i]) .& (df_full.θ .== θ[j]),:]
    df_single = df_single[abs.(df_single.E) .< 0.3,:]
    val, idx = findmax(df_single.PN_mean)
    df_temp = DataFrame(E = df_single[idx, :E], V = df_single[idx,:V], θ = df_single[idx,:θ], PN_mean = val, PN_err = df_single[idx,:PN_err],N_E = df_single[idx,:N_E])
    append!(dfmax,  df_temp)
end

dfmax_del0 = DataFrame(E = Float64[], V = Float64[], θ = Float64[], PN_mean = Float64[], PN_err = Float64[], N_E = Int64[])
for j in 1:length(θ)
    df_single = df_del0[df_del0.θ .== θ[j],:]
    df_single = df_single[abs.(df_single.E) .< 0.3,:]
    val, idx = findmax(df_single.PN_mean)
    df_temp = DataFrame(E = df_single[idx, :E], V = -3, θ = df_single[idx,:θ], PN_mean = val, PN_err = df_single[idx,:PN_err],N_E = df_single[idx,:N_E])
    append!(dfmax_del0,  df_temp)
end

dfmax_delinf = DataFrame(E = Float64[], V = Float64[], θ = Float64[], PN_mean = Float64[], PN_err = Float64[], N_E = Int64[])
for j in 1:length(θ)
    df_single = df_delinf[df_delinf.θ .== θ[j],:]
    df_single = df_single[abs.(df_single.E) .< 0.3,:]
    val, idx = findmax(df_single.PN_mean)
    df_temp = DataFrame(E = df_single[idx, :E], V = 6, θ = df_single[idx,:θ], PN_mean = val, PN_err = df_single[idx,:PN_err],N_E = df_single[idx,:N_E])
    append!(dfmax_delinf,  df_temp)
end

df_fig2_1 = dfmax[(dfmax.θ .== 0.25), :]
df_fig2_2 = dfmax[(dfmax.V .== dfmax.V[100]) .| (dfmax.V .== dfmax.V[200]) .| (dfmax.V .== dfmax.V[400]) .| (dfmax.V .== dfmax.V[600]), :]
df_fig2_2.V = round.(df_fig2_2.V, digits = 0)

cscheme = ColorSchemes.jet1

## PN at a fixed E
Gadfly.with_theme(theme) do
    set_default_plot_size(13cm, 10cm)
#------------------------------ PN vs θ ------------------------------#
    xticks = collect(0.00:0.05:0.25)
    yticks = 1:2:10
    l1 = layer(df_fig2_2, x = :θ, y = :PN_mean, color = :V, Geom.line, Geom.point, order = 1)
    l2 = layer(dfmax_del0, x = :θ, y = :PN_mean, color =[colorant"black"], Geom.line, order = 2, style(line_style = [:dash], line_width = 2.5pt))
    l3 = layer(dfmax_delinf, x = :θ, y = :PN_mean, color = [colorant"black"], Geom.line, order = 3,  style(line_style = [:dash], line_width = 2.5pt))
    p1 = plot(l1, l2, l3,
        Guide.xlabel("Θ/π"),
        Guide.ylabel("⟨PN⟩"),
        Guide.xticks(ticks = xticks),
        Guide.yticks(ticks = yticks),
        Guide.colorkey(title="Log<sub>10</sub>[δ]", pos=[0.22w,-0.24h]),
        Guide.annotation(compose(context(), text(0.55w, 0.2h, "δ = ∞"),fontsize(20pt))),
        Guide.annotation(compose(context(), text(0.63w, 0.95h, "δ = 0"),fontsize(20pt))),
        color = ["black"],
        style(point_size = 3.5pt, line_width = 2pt),
    );
    display(p1)
    draw(PDF(savedir*"p1_fig.pdf"), p1)

end


#------------------------------ PN vs V ------------------------------#
Gadfly.with_theme(theme) do
    set_default_plot_size(13cm, 10cm)
    xticks = -2:1:5
    yticks = 2:2:10
    l_del0 = layer(dfmax_del0[dfmax_del0.θ .== 0.25, :], x = :V, y = :PN_mean, color =[colorant"black"], Geom.point, order = 2)
    l_delinf = layer(dfmax_delinf[dfmax_delinf.θ .== 0.25, :], x = :V, y = :PN_mean, color =[colorant"black"], Geom.point, order = 3)

    l2 = layer(df_fig2_1, x =:V, y =:PN_mean, Geom.line, Geom.point, color = [colorant"black"])
    p2 = plot(l_del0,l_delinf, l2,
        Guide.xlabel("Log<sub>10</sub>[δ]"),
        Guide.ylabel("⟨PN⟩"),
        Guide.xticks(ticks = xticks),
        Guide.yticks(ticks = yticks),
        Guide.annotation(compose(context(), text(0.5w, 0.2h, "Θ = 0.25π"),fontsize(22pt))),
        style(line_width = 2pt, point_size = 3.5pt)

    );
    display(p2)
    draw(PDF(savedir*"p2_fig.pdf"), p2)

end

A = copy(dfmax.PN_mean)
A = collect(reshape(A, 25,25))
xs = range(p.V_min, p.V_max, length = p.V_num)
ys = range(p.θ_min, p.θ_max, length = p.θ_num)

x_interp = range(xs[1], xs[end], length = 100)
y_interp = range(ys[1], ys[end], length = 100)
interp_linear = LinearInterpolation((xs, ys), A)
PN_interp = interp_linear(x_interp, y_interp)
PN_interp = reshape(PN_interp, 100*100)
xy_interp = reshape(collect(Iterators.product(y_interp,x_interp)), 100*100)
V_interp = getindex.(xy_interp, 2)
θ_interp = getindex.(xy_interp ,1)


df_interp =DataFrame(V = V_interp, θ = θ_interp, PN_mean = PN_interp)
Gadfly.with_theme(theme) do
    set_default_plot_size(15cm, 10cm)
#------------------------------ PN vs θ, DELTA ------------------------------#
    yticks = collect(0.00:0.05:0.25)
    yticks[1] = 0.01
    xticks = -2:1:5
    p3 = plot(df_interp, color = :PN_mean, y = :θ, x = :V, Geom.rectbin,
        Scale.ContinuousColorScale(p -> get(ColorSchemes.viridis, p)),
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.ylabel("Θ/π"),
        Guide.xlabel("Log<sub>10</sub>[δ]"),
        Guide.colorkey(title="⟨PN⟩"),
        )
    display(p3)
    draw(PDF(savedir*"p3_fig.pdf"), p3)

end

using JSON
using ColorSchemes
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
binnum = 101
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
    df_single = df_single[abs.(df_single.E) .< 0.2,:]
    val, idx = findmax(df_single.PN_mean)
    df_temp = DataFrame(E = df_single[idx, :E], V = df_single[idx,:V], θ = df_single[idx,:θ], PN_mean = val, PN_err = df_single[idx,:PN_err],N_E = df_single[idx,:N_E])
    append!(dfmax,  df_temp)
end
## PN at a fixed E
Gadfly.with_theme(theme) do
    cscheme = ColorSchemes.jet1
#------------------------------ PN vs θ ------------------------------#
    xticks = collect(0.00:0.05:0.25)
    yticks = 1:2:9
    l1 = layer(dfmax, x = :θ, y = :PN_mean, color=:V, Geom.line, order = 1)
    p1 = plot(l1,
        Guide.xlabel("Θ/π"),
        Guide.ylabel("PN"),
        Guide.xticks(ticks = xticks),
        Guide.yticks(ticks = yticks),
        Guide.colorkey(title="log<sub>10</sub>[δ]")
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
    );
    display(p1)
#------------------------------ PN vs V ------------------------------#
    xticks = -2:1:5
    yticks = 0:2.5:10
    p2 = plot(dfmax, x =:V, y =:PN_mean, color=:θ, Geom.line,
        Guide.xlabel("log<sub>10</sub>[δ]"),
        Guide.ylabel("PN"),
        Guide.xticks(ticks = xticks),
        Guide.yticks(ticks = yticks),
        Guide.colorkey(title="Θ/π"),
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
    );
    display(p2)
#------------------------------ PN vs θ, DELTA ------------------------------#
    xticks = collect(0.00:0.05:0.25)
    xticks[1] = 0.01
    yticks = -2:1:5
    p3 = plot(dfmax, color = :PN_mean, y = :V, x = :θ, Geom.rectbin,
        Scale.ContinuousColorScale(p -> get(ColorSchemes.viridis, p)),
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("Θ/π"),
        Guide.ylabel("log<sub>10</sub>[δ]"),
        Guide.colorkey(title="PN"),
        )
    display(p3)
#------------------------------ PN vs E at Θ=0.25π ------------------------------#
    xticks = -0.5:0.2:0.5
    yticks = 0:2.5:10
    p4 = plot(dfth_full, x =:E, y =:PN_mean, color=:V, Geom.line,
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("E"),
        Guide.ylabel("PN"),
        Guide.colorkey(title="log<sub>10</sub>[δ]"),
        style(line_width=0.5mm)
        )
    display(p4)

#------------------------------ PN vs DELTA at Θ=0.25π------------------------------#
    xticks = -2:1:5
    yticks = 0:2.5:10
    p5 = plot(dfth_full[-0.4 .< dfth_full.E .< 0.0, :], x =:V, y =:PN_mean, color=:E, Geom.line,
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("log<sub>10</sub>[δ]"),
        Guide.ylabel("PN"),
        Guide.colorkey(title="E"),
        style(line_width=0.5mm)
        )
    display(p5)
#------------------------------ PN vs E, DELTA at Θ=0.25π ------------------------------#
    xticks = -0.5:0.2:0.5
    yticks = -2:1:5
    p6 = plot(dfth_full, x =:E, y =:V, color=:PN_mean, Geom.rectbin,
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("E"),
        Guide.ylabel("log<sub>10</sub>[δ]"),
        Guide.colorkey(title="PN")
        )
    display(p6)

    draw(PDF(savedir*"p1.pdf"), p1)
    draw(PDF(savedir*"p2.pdf"), p2)
    draw(PDF(savedir*"p3.pdf"), p3)
    draw(PDF(savedir*"p4.pdf"), p4)
    draw(PDF(savedir*"p5.pdf"), p5)
    draw(PDF(savedir*"p6.pdf"), p6)
end

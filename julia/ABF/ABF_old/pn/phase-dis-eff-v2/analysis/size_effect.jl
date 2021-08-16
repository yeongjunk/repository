using JSON
using ColorSchemes
using Compose
include("../library/abf2_analysis_tool.jl")
include("../library/abf2_pnscan.jl")

fullscan_dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-nu2-sf-size-effect/rawdata"
arr_fn = readdir_jld(fullscan_dir, join = true)
arr_fn_short = readdir_jld(fullscan_dir, join = false)

save_dir = "/Users/pcs/data/ABF/analysis/nu2-pn-phase-dis-eff-new/size-effect/"

binnum = [3001;1001;301;101]
E_edges = [range(-0.5, 0.5, length = binnum[i]+1) for i in 1:length(binnum)]

## Vary the size of the uniform energy bins
## Apr4 stopped at i = 182
set_default_plot_size(26cm, 20cm)
Gadfly.push_theme(theme)

L = [1001; 101; 3001; 301]
## Find maximum
dfmax = DataFrame(E = Float64[], PN_max = Float64[],
    PN_err = Float64[], PN_std = Float64[], L = Int64[], binnum = Int64[])
for j in 1:length(E_edges)
    for i in 1:4
        df = summary_single(arr_fn[i], E_edges[j])
        df = df[abs.(df.E) .< 0.05, :]
        val, idx = findmax(df.PN_mean)
        push!(dfmax, (df[idx,:E], df[idx,:PN_mean], df[idx,:PN_err],df[idx,:PN_std], L[i], binnum[j]))
    end
    println(j)
end
dfmax.L_inv = 1 ./dfmax.L
dfmax.errmax = dfmax.PN_max .+ dfmax.PN_err
dfmax.errmin = dfmax.PN_max .- dfmax.PN_err
sort!(dfmax, :L, );

p = plot(dfmax, x = :L_inv, y = :PN_max,
    ymax = :errmax, ymin =:errmin, color = :binnum,
    Geom.line, Geom.point, Geom.errorbar,
    Scale.color_discrete_hue(), style(point_size = 3pt),
    Guide.xlabel("1/L"),
    Guide.ylabel("max(⟨PN⟩)"))

draw(PNG(save_dir*"result1.png", dpi = 150), p)

p2 = plot(dfmax, x = :binnum, y = :PN_max, color = :L,
    Geom.line,Geom.point, Scale.color_discrete_hue())

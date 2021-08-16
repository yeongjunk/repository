using JSON
using ColorSchemes
using Compose
include("../library/abf2_analysis_tool.jl")
include("../library/abf2_pnscan.jl")

fullscan_dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-nu2-sf-size-effect/rawdata"
#fullscan_dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-abf-sf/merged"
arr_fn = readdir_jld(fullscan_dir, join = true)
arr_fn_short = readdir_jld(fullscan_dir, join = false)

save_dir = "/Users/pcs/data/ABF/analysis/nu2-pn-phase-dis-eff-new/binsize-effect/"

binnums = [3001;1001;301;101; 31]
E_edges = [range(-0.5, 0.5, length = i+1) for i in binnums]


## Vary the size of the uniform energy bins
## Apr4 stopped at i = 182
set_default_plot_size(26cm, 20cm)
Gadfly.push_theme(theme)
for i in 1:4
    dict = JLD.load(arr_fn[i])
    θ = dict["θ"]
    V = dict["V"]
    df_raw = dict["data"]
    V_log = round(log10(V), digits = 1)

    df_raw = df_raw[(abs.(dict["data"].E) .< 0.01), :]
    df_mid = midband(dict["data"], 10^-6)
    pn_mid = mean(df_mid.PN)
    pn_mid_err = std(df_mid.PN)/sqrt(size(df_mid,1))
    println("band center finished")

    df = DataFrame(E = Float64[], V = Float64[], θ=Float64[], PN_mean = Float64[], PN_err = Float64[], PN_std = Float64[], N_E = Float64[], binnum = String[], Bandwidth = Float64[])
    for j in 1:length(E_edges)
        print("Binning: $(j)/$(length(E_edges))\r")
        df_temp = summary_single(V, θ, df_raw.E, df_raw.PN, E_edges[j])
        newrow = Array{String}(undef, size(df_temp,1))
        fill!(newrow, "$(length(E_edges[j])-1)")
        df_temp.binnum = newrow
        append!(df, df_temp)
    end

    df_mid_trend = DataFrame(binnum = Int64[], PN_mean = Float64[], PN_max = Float64[], PN_min = Float64[], PN_peak = Float64[], PN_max2 = Float64[], PN_min2 = Float64[],)

    for j in 1:length(binnums)
        df_temp = df[(df.binnum .== "$(binnums[j])") .& (abs.(df.E) .< 0.3), :]
        peakval,idx = findmax(df_temp.PN_mean[.!isnan.(df_temp.PN_mean)])

        z_PN_mean = df_temp.PN_mean[binnums[j]÷2+1]
        z_PN_max = z_PN_mean + df_temp.PN_err[binnums[j]÷2+1]
        z_PN_min = z_PN_mean - df_temp.PN_err[binnums[j]÷2+1]

        PN_peak = peakval
        PN_max2 = peakval + df_temp.PN_err[idx]
        PN_min2 = peakval - df_temp.PN_err[idx]

        append!(df_mid_trend, DataFrame(binnum = [binnums[j]], PN_mean = [z_PN_mean], PN_max = [z_PN_max], PN_min = [z_PN_min],
            PN_peak = [PN_peak], PN_max2 = [PN_max2], PN_min2 = [PN_min2]))
    end
    println("Plotting")
    println(maximum(df.E))
    p1 = plot(df, x = :E, y = :PN_mean, color = :binnum, Geom.line,
        style(line_width = 1.2pt,), Guide.ylabel("⟨PN⟩"))

    p1_inset = plot(df[abs.(df.E) .< 0.01, :], x = :E, y = :PN_mean, color = :binnum, Geom.line,
        style(line_width = 1.2pt, key_position = :none, panel_fill=colorant"CornSilk"), Guide.ylabel("⟨PN⟩"));

    p1_inset2 = plot(df_mid_trend, x = :binnum, y = :PN_peak, ymax = :PN_max2, ymin = :PN_min2, Geom.line, Geom.point, Geom.errorbar, Scale.x_log10,
        style(highlight_width=0.0pt, key_position = :none, panel_fill=colorant"CornSilk"),
        Guide.ylabel("max(⟨PN⟩)"),
        Guide.xticks(ticks = 1:1:4))

    l1 = layer(x = [0.], y = [pn_mid], color = ["near E=0"], ymax = [pn_mid + pn_mid_err], ymin = [pn_mid - pn_mid_err],
        Geom.point,
        Geom.errorbar,
        style(highlight_width=0.0pt, point_size = 3pt)
        )

    p_full = plot();
    push!(p1_inset, l1);
    p_full = push!(p1, l1,
        Guide.annotation(compose(context(), text(0.1w, 0.1h, "log<sub>10</sub>[δ] = $(V_log) , θ = $(θ)π"))),
        Guide.annotation(compose(context(.55w, 0.02h ,0.45w, 0.45h), render(p1_inset))),
        Guide.annotation(compose(context(.55w, 0.55h ,0.45w, 0.45h), render(p1_inset2))),
        )
    ## Non-uniform energy bins near zero
    draw(PNG(save_dir*"$(arr_fn_short[i]).png", dpi = 150), p_full)
    println("$(i)/$(length(arr_fn))")
end
Gadfly.pop_theme()

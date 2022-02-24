using JLD
using DataFrames
using Statistics
using Gadfly
using StatsBase
import Cairo, Fontconfig
using Compose
using LinearAlgebra


## Reading Files
@doc """
read directory, and find all files that include ".jld"
"""
function readdir_jld(dir; join = false, sort = false)
    s = filter!(s->occursin(r".jld", s), readdir(dir, join = join, sort = sort))
end

## Process data
@doc """
Find zero energy states, |E| < tol
"""
function midband(df, tol)
    return midband_inner(df.E, df.PN, tol)
end

@doc """
Barriar function for midband
"""
function midband_inner(E, PN, tol)
    idx = findall(x-> abs(x) < tol, E)
    return DataFrame(E = E[idx], PN = PN[idx])
end

@doc """
Compute the bandwidth
"""
function bandwidth(E::AbstractArray{Float64})
    return Float64(maximum(E) - minimum(E))
end

@doc """
Normalize bandwidth
"""
function normalize_bw!(E::AbstractArray{Float64}; bw = nothing)
    if isnothing(bw)
        bw = bandwidth(E)
    end
        E .= E/bw
end

@doc """
binary search binning x (non-uniform edges, unsorted array)
"""
function bin_x(x,x_edges)
    L = 1 #leftmost edge
    R = length(x_edges) #rightmost edge
    #j = 1;
    while L+1 < R
        m = (L+R)÷2
        # println("$(j),$(L), $(m), $(R)")
        if x  >= x_edges[m]
             L = m
        elseif x < x_edges[m]
            R = m
        else
            return 0
        end
        #j += 1
    end
    return L
end

function bin_arr_sorted(x_arr, x_edges)
    i = 2
    j = 1
    lbl = zeros(Int64, length(x_arr))
    for i in 2:length(x_edges)
        while (j <= length(x_arr)) && (x_edges[i-1] <= x_arr[j] < x_edges[i])
            lbl[j] = i-1
            j += 1
        end
    end
    return lbl
end


@doc """
label array of data by bin index
"""
function bin_arr(x_arr, x_edges)
    lbl = bin_x.(x_arr, Ref(x_edges))
    return lbl
end

@doc """
Summarize a single scan result
"""
function summary_single(V, θ, E, PN, E_edges; ztol = nothing)
    bw = bandwidth(E)
    normalize_bw!(E, bw=bw)

    binnum = bin_arr(E, E_edges)

    #Statistics
    d_avg = Float64[];
    d_err = Float64[];
    d_std = Float64[];
    d_n = Int64[];

    for i in 1:(length(E_edges)-1)
        PN_at_bin_i = PN[binnum .== i]
        PN_std = std(PN_at_bin_i)
        n = size(PN_at_bin_i,1)
        push!(d_n, n)
        push!(d_avg, mean(PN_at_bin_i))
        push!(d_err, PN_std/sqrt(n))
        push!(d_std, PN_std)
    end

    df_sum = DataFrame(E = midpoints(E_edges), V = V, θ = θ, PN_mean = d_avg, PN_std = d_std, PN_err = d_err, N_E =d_n, Bandwidth = bw)

    return df_sum
end


@doc """
Summarize a single scan result
"""
@inline function summary_single(fn::String, E_edges; ztol = nothing)
    data = JLD.load(fn)
    df = data["data"]
    df_sum = summary_single(data["V"], data["θ"], data["data"].E, data["data"].PN, E_edges; ztol)
    return df_sum
end


## Gadfly plot theme
ft = "CMU Serif";
theme = Theme(
    panel_fill= "white",
    #-----FONT TYPE-----#
    major_label_font= ft,
    minor_label_font= ft,
    key_title_font  = ft,
    key_label_font  = ft,
    point_label_font = ft,
    #-----FONT SIZE-----#
    major_label_font_size=24pt,
    minor_label_font_size=18pt,
    key_title_font_size=18pt,
    key_label_font_size=18pt,
    point_label_font_size = 18pt,
    #-----FONT COLOR-----#
    minor_label_color = "black",
    major_label_color = "black",
    key_title_color  =  "black",
    key_label_color  =  "black",
    #-----PANELS & GRID----#
    panel_stroke    =   "black",
    panel_line_width = 2pt,
    grid_line_width = 1pt,
    grid_line_style = :solid,
    plot_padding = [1mm, 5mm, 1mm, 1mm],
    line_width=0.5mm
)

@doc """
Summarize a single scan result with plots
"""
function summary_single_visual(fn, E_edges, theme)
    data = JLD.load(fn)
    V = data["V"]
    θ = data["θ"]
    V_log = round(log10(data["V"]), digits = 2)
    df = data["data"]
    bw = bandwidth(df.E)
    normalize_bw!(df.E, bw = bw)
    df.binnum = bin_arr(df.E, E_edges)

    #Statistics
    d_avg = Float64[];
    d_err = Float64[];
    d_std = Float64[];
    d_n = Int64[];

    for i in 1:(length(E_edges)-1)
        PN_at_bin_i = df[df.binnum .== i,:].PN
        PN_std = std(PN_at_bin_i)
        n = size(PN_at_bin_i,1)
        push!(d_n, n)
        push!(d_avg, mean(PN_at_bin_i))
        push!(d_err, PN_std/sqrt(n))
        push!(d_std, PN_std)
    end

    df_sum = DataFrame(E = midpoints(E_edges), V = V, θ = θ, PN_mean = d_avg, PN_std = d_std, PN_err = d_err, N_E = normalize(d_n), Bandwidth = bw)



    xticks = -0.5:0.25:0.5
    xticks2 = -0.005:0.005:0.005
    Gadfly.with_theme(theme) do
        dfsp = df[1:50000,:] #sparse
        p1 = plot(dfsp, x = :E, y = :PN, color = :binnum,
        Scale.color_discrete, style(key_position = :none, highlight_width=0.5pt, point_size = 1.3pt),
        Guide.annotation(compose(context(), text(0.1w, 0.25h, "log<sub>10</sub>[δ] = $(V_log) , θ = $(θ)π"))),
        Guide.xticks(ticks = xticks))

        l1 = layer(df_sum, x = :E, y = :PN_mean, ymax = df_sum.PN_mean .+ df_sum.PN_err, ymin = df_sum.PN_mean .- df_sum.PN_err,
        Geom.line, Geom.point, Geom.errorbar)
        df_mid = midband(df, 10^-5)
        l2 = layer(x = [0], y = [mean(df_mid.PN)], ymax =[mean(df_mid.PN)] + [std(df_mid.PN)/sqrt(size(df_mid,1))], ymin = [mean(df_mid.PN)] - [std(df_mid.PN)/sqrt(size(df_mid,1))], Geom.point, Geom.errorbar, color = [colorant"red"])

        p2 = plot(l1,l2, style(highlight_width=0pt),
        Guide.xticks(ticks = xticks),
        Guide.xlabel("E"),Guide.ylabel("mean(PN)"),
        style(key_position = :none, highlight_width=0pt, point_size = 1.5pt))

        p3 = plot(dfsp[abs.(dfsp.E) .< 0.005, :], x = :E, y = :PN, color = :binnum,
         style(key_position = :none, highlight_width=0.5pt, point_size = 1.5pt),
        Guide.xticks(ticks = xticks2))

        p4 = plot(df_sum, x = :E, y = :N_E, Geom.line, Geom.point,
        Guide.xticks(ticks = xticks),
        Guide.xlabel("E"),Guide.ylabel("P(E)"),  style(key_position = :none, highlight_width=0.5pt, point_size = 1.5pt))

        set_default_plot_size(24cm, 16cm)
        p_full = gridstack([p1 p2; p3 p4])
        return p_full
    end
end


function summary_full(arr_fn::AbstractArray{String}, E_edges; exc_zero = false)
     df_full = DataFrame(E = Float64[], V = Float64[], θ = Float64[], PN_mean = Float64[], PN_err = Float64[], N_E = Int64[])
     for i in 1:length(arr_fn)
         print("$(i)/$(length(arr_fn))\r")
         df_single = summary_single(arr_fn[i], E_edges; exc_zero)
         append!(df_full, df_single)
     end
     return df_full
 end

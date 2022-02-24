using DataFrames
using JLD
using StatsBase
using CSV
include("../library/abf2_pnscan.jl")
topdir = "/Users/pcs/data/ABF/rawdata/nu2-pn/phase-dis-eff2/"
fn = "pn_merged.jld"


function reshape1d(data)
    return reshape(data, length(data))
end
function load_data(data::Dict)
    return data["PN_mean"], data["PN_var"], data["num_samp"], data["PN_hist"],data["E_bw"]
end
function coord_column(x,y)
    domain = collect(Iterators.product(x,y))
    domain_col = reshape(domain, length(domain))
    x_col = getindex.(domain_col,1)
    y_col = getindex.(domain_col,2)
    return hcat(x_col,y_col)
end
function coord_column(x,y,z)
    domain = collect(Iterators.product(x,y,z))
    domain_col = reshape(domain, length(domain))
    x_col = getindex.(domain_col,1)
    y_col = getindex.(domain_col,2)
    z_col = getindex.(domain_col,3)
    return hcat(x_col,y_col, z_col)
end
function dframe_pn(data::Dict)
    p = readconfig(data)
    θ, δ, E = expand_params(p)
    pn_mean, pn_var, N, pn_hist, E_bw = load_data(data)

    domain_col = coord_column(θ, δ, E)
    if p.δ_log == true
        domain_col[:,2] = log10.(domain_col[:,2])
    end

    PN_mean_1d = reshape1d(pn_mean)
    PN_var_1d = reshape1d(pn_var)
    N_1d = reshape1d(N)

    if p.δ_log == true
        df = DataFrame(
            Theta = domain_col[:,1], L10_Delta = domain_col[:,2],E = domain_col[:,3],
            PN_mean = PN_mean_1d,
            PN_var = PN_var_1d,
            N_sample = N_1d,
            )
    else
        df = DataFrame(
            Theta = domain_col[:,1], Delta = domain_col[:,2],E = domain_col[:,3],
            PN_mean = PN_mean_1d,
            PN_var = PN_var_1d,
            N_sample = N_1d,
        )
    end
    return df
end
function dframe_bw(data::Dict)
    p = readconfig(data)
    θ, δ, _ = expand_params(p)
    _, __, ___, ____, E_bw = load_data(data)

    domain = coord_column(θ,δ)
    if p.δ_log == true
        domain[:,2] = log10.(domain[:,2])
    end
    E_bw_1d = reshape1d(E_bw)
    if p.δ_log == true
        df = DataFrame(
            Theta = domain[:,1], L10_Delta = domain[:,2],
            Δ = E_bw_1d,
        )
    else
        df = DataFrame(
            Theta = domain[:,1], Delta = domain[:,2],
            Δ = E_bw_1d,
        )
    end
end
#-------------------------LOAD-------------------------#
data = JLD.load(topdir*fn)
data["E_abs_bw"] = false
df = dframe_pn(data)
df2 = dframe_bw(data)
#-------------------------SAVE DF-------------------------#
CSV.write(topdir*"pn_dataframe.csv", df)
CSV.write(topdir*"pn_bw.csv", df2)
#-------------------------HISTOGRAM-------------------------#
p = readconfig(data)
θ, δ, _ = expand_params(p)
PN_hist = data["PN_hist"]
#Extremely unneccessary histogram bulk
for i in 1:length(θ), j in 1:length(δ)
    df_hist_arr = [0. 0. 0. 0. 0.];
    xaxis = round.(collect(midpoints(pn_hist[i,j].edges[1])), digits = 12)
    yaxis = midpoints(pn_hist[i,j].edges[2])
    xyaxis = Iterators.product(xaxis,yaxis)
    xyaxis = reshape(collect(xyaxis), length(xyaxis))
    idx_1 = repeat([θ[i]], outer = length(xyaxis))
    idx_2 = repeat([δ[j]], outer = length(xyaxis))
    idx_3 = getindex.(xyaxis, 1)
    idx_4 = getindex.(xyaxis, 2)
    val = reshape(pn_hist[i,j].weights, length(xyaxis))
    df_hist_arr = vcat(df_hist_arr, hcat(idx_1, idx_2, idx_3, idx_4, val))
    df_hist_arr = Array{Float64}(df_hist_arr[2:end,:])

    df_hist = DataFrame(
        Theta = df_hist_arr[:,1],
        L10_Delta = df_hist_arr[:,2],
        E = df_hist_arr[:,3],
        PN = df_hist_arr[:,4],
        Num = convert.(Int64, df_hist_arr[:,5])
    )
    CSV.write(topdir*"hist-csv/"*"th$(i)-del$(j).csv", df_hist)
    println("saved-th$(i)-del$(j)")
end

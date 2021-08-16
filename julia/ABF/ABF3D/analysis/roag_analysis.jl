using DataFrames, CSV
const modulepath =  "/Users/pcs/codes/chain/ABF/module"
push!(LOAD_PATH, modulepath)
using ROAG, Binning
using Statistics
using StatsBase
using MAT
## Function definitions
function load_file(dir, i, j, L, R, θ)
    subdir = "L$(L[i])/"
    fn = "L$(L[i])_Th$(j)_R$(R[i]).csv"
    return Tables.columntable(CSV.read(dir*subdir*fn, DataFrame))
end

function load_file(dir, subsubdir, i, j, L, R, θ)
    subdir = "L$(L[i])/"
    fn = "L$(L[i])_Th$(j)_R$(R[i]).csv"
    return Tables.columntable(CSV.read(dir*subdir*subsubdir*fn, DataFrame))
end

function roag_mean_binned(E::AbstractArray{T}, E_edges::AbstractArray{T}) where T <: Number
    r_mean = Vector{Float64}(undef, length(E_edges)-1)
    lbl = binning_unsorted(E, E_edges) # binning
    for i in 1:(length(E_edges) - 1)
        r = E[findall(a -> a == i, lbl)]
        if length(r) < 3
            r_mean[i] = NaN
        else
            roag!(r)
            r_mean[i] = mean(r)
        end
    end
    return r_mean::Vector{Float64}
end

function roag_scan(;dir, L, R, θ, E_edges)
    r_m_m = Array{Float64}(undef, length(E_edges)-1 ,length(θ), length(L))
    r_m_std = similar(r_m_m)
    for i in 1:length(L)
        println("i = ", i)
        for j in 1:length(θ)
            df = load_file(dir, i, j, L, R, θ)
            r_m = Array{Float64}(undef, length(E_edges)-1, R[i])
            for r in 1:R[i]
                idx = searchsorted(df.r, r)
                E = df.E[idx]
                r_m[:, r] .= roag_mean_binned(E, E_edges)
            end
            for k in 1:length(E_edges) - 1
                r_m_m[k, j, i] = mean(filter(!isnan, r_m[k, :]))
                r_m_std[k, j, i] = std(filter(!isnan, r_m[k, :]))
            end
        end
    end
    return r_m_m, r_m_std
end


function roag_scan(;dir, L, R, θ, E_edges, subsubdir = "")
    r_m_m = Array{Float64}(undef, length(E_edges)-1 ,length(θ), length(L))
    r_m_std = similar(r_m_m)
    for i in 1:length(L)
        println("i = ", i)
        for j in 1:length(θ)
            df = load_file(dir, subsubdir, i, j, L, R, θ)
            r_m = Array{Float64}(undef, length(E_edges)-1, R[i])
            for r in 1:R[i]
                idx = searchsorted(df.r, r)
                E = df.E[idx]
                r_m[:, r] .= roag_mean_binned(E, E_edges)
            end
            for k in 1:length(E_edges) - 1
                r_m_m[k, j, i] = mean(filter(!isnan, r_m[k, :]))
                r_m_std[k, j, i] = std(filter(!isnan, r_m[k, :]))/sqrt(length(filter(!isnan, r_m[k, :])))
            end
        end
    end
    return r_m_m, r_m_std
end

## Parameters
L = [10 15 20 25]
R = [6000 1800 600 300]
θ = range(0.0001, 0.25, length = 26)
E_edges = collect(0.65:0.1:1.35)
pushfirst!(E_edges,0.625)
pushfirst!(E_edges, 0.6)
push!(E_edges,  1.375)
push!(E_edges,  1.4)
push!(E_edges, 1.30)
push!(E_edges, 0.70)
sort!(E_edges)
## Parameters/fine
L_f = [10 15 20 25]
R_f = [10000 4000 1600 600]
θ_f = range(0.09, 0.11, length = 11)
E_edges_f = collect(range(0.95, 1.05, length = 2))

##
dir = "/Users/pcs/data/ABF3D/full-e/"
dir_f = "/Users/pcs/data/ABF3D/full-e-fine/"
savedir = "/Users/pcs/data/ABF3D/analyzed-data/"
##
r_m_m, r_m_std = roag_scan(dir = dir, subsubdir = "", L = L, R = R , θ = θ, E_edges = E_edges)
r_m_m_f, r_m_std_f = roag_scan(dir = dir_f,subsubdir = "", L = L_f, R = R_f , θ = θ_f, E_edges = E_edges_f)


r_m_m_c = r_m_m[7,:,1:4]
idx = [collect(1:9); collect(13:26)]
θ_c = θ[idx]
r_m_m_c = r_m_m_c[idx, :]
θ_merged = vcat(θ_c[1:9], collect(θ_f), θ_c[10:end])
r_m_m_merged = Array{Float64}(undef, length(θ_merged), 4)
r_m_m_merged[1:9, :] = r_m_m_c[1:9, :]
r_m_m_merged[10:20, :] = r_m_m_f
r_m_m_merged[21:end, :] = r_m_m_c[10:end,:]


matwrite(savedir*"roag_merged.mat",
    Dict("r_m_m" => r_m_m,
    "r_m_std" => r_m_std,
    "r_m_m_f" => r_m_m_f,
    "r_m_std_f" => r_m_std_f,
    "th_merged" => θ_merged,
    "r_m_m_merged" => r_m_m_merged)) # note that this does NOT introduce a variable ``varname`` into scope

r_m_data = matread(r_m_dir)
r_m_m_f = r_m_data["r_m_m_f"]
r_m_std_f = r_m_data["r_m_std_f"]

df = [DataFrame(th = θ_f, r = r_m_m_f[1,:,i], std = r_m_std_f[1,:,i]) for i in 1:length(L)]
for i in 1:length(L)
    CSV.write(savedir*"ABF3D_sf_r_L$(L[i]).csv", df[i])
end

using MAT
using LinearAlgebra
using Plots
using DataFrames
using CSV
@doc """
A struct where parameters are stored. This will be used as inputs.
"""
struct Parameters
    E_min::Float64
    E_max::Float64
    E_num::Int64

    W_min::Float64
    W_max::Float64
    W_num::Int64

    θI_min::Float64
    θI_max::Float64
    θI_num::Int64

    θII_min::Float64
    θII_max::Float64
    θII_num::Int64

    seed::Int64
    N::Int64
end

@doc """
Read configuration JSON file, and construct a Parameters struct from it.
"""
function config_read(config::Dict{String,Any})
    p = Parameters(config["E_min"],
           config["E_max"],
           config["E_num"],

           config["W_min"],
           config["W_max"],
           config["W_num"],

           config["theta_min"][1],
           config["theta_max"][1],
           config["theta_num"][1],

           config["theta_min"][2],
           config["theta_max"][2],
           config["theta_num"][2],

           config["seed"],
           config["N"],
           )

           return p
end

function params_expand(p::Parameters)
    E_all   = collect(LinRange(p.E_min, p.E_max, p.E_num))     #   Range of eigenenergies
    W_all   = collect(LinRange(p.W_min, p.W_max, p.W_num))     #   Range of disorder strengths

    #   Ranges of angles parameterising the nu=2 ABF models
    θI_all  = collect(LinRange(p.θI_min, p.θI_max, p.θI_num));
    θII_all  = collect(LinRange(p.θII_min, p.θII_max, p.θI_num));

    # r_all   = rand_array_create(p.seed, p.N)              # pregenerate random numbers

    return E_all, W_all, θI_all, θII_all
end

@doc """
x,y axis are the first and the second axis of matrix z, respectively
"""
function data_df_style(x::AbstractArray{Float64, 1}, y::AbstractArray{Float64, 1}, z::AbstractArray{Float64, 2})
    z_reshape = reshape(z, size(z, 1)*size(z, 1), )
    x_reshape = repeat(x, outer = length(y))
    y_reshape = repeat(y, inner = length(x))
    return x_reshape, y_reshape, z_reshape
end

@doc """
read ξ and parameters from saved .mat file
"""
function ξ_read(matfilename::String)
    vars = matread(matfilename) #Dictionary

    # Localization length
    ξ_all   = vars["xi_all"]
    params = config_read(vars) #create a Parameters struct

    return ξ_all, params
end

xi_all, params = ξ_read("/Users/pcs/data/ABF-sum/1d-full-on-tm-xi.mat")
E, W, th1, th2 = params_expand(params)
xi_diag = Array{Float64}(undef, 51, 51, 51)
for i in 1:51
    xi_diag[:, :, i] = xi_all[:, :, i, i]
end

th_val = [pi/4, pi/8, pi/16, 0]
th_val_s = ["0.25", "0.125", "0.0625", "0"]
for i in 1:length(th_val)
    val, idx = findmin(abs.(th1 .- th_val[i]))
    E_rs, W_rs, xi_rs = data_df_style(E, W, xi_diag[:, :, idx])
    df = DataFrame(E = round.(E_rs, sigdigits = 12), W = W_rs, xi = round.(xi_rs, sigdigits = 12))
    CSV.write("/Users/pcs/data/ABF-sum/finalized-data-csv/xi-$(th_val_s[i]).csv", df)
end

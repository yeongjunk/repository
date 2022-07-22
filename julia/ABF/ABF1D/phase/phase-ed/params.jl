using Parameters

@with_kw struct Params
    l::Vector{Int64} = [1]
    L::Int64 = 100
    θ::Vector{Float64} = [0.25]
    W::Vector{Float64} = [0.]
    E::Vector{Float64} = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0] # Normalized energy, [0, 1]]
    q::Vector{Float64} = [2, 3]
    seed::Int64 = 1234
    R::Int64 = 100
    nev::Int64  = 10
    V::Float64 = 1.
    E_window_width::Int64 = 200 # The full width will be 4/E_window_width, of energy band [-1, 1]
    num_blas::Int64 = 1 
    energy_path::String
end

function readconfig(config::Dict)
    if typeof(config["l"]) <: Integer
        l = [config["l"]];
    else
        l = config["l"]
    end
    L = config["L"]
    θ = Vector(config["th"]) 
    W = Vector(config["W"])
    E = Vector(config["E"])
    q = isa(config["q"], AbstractArray) ? Array{Int64}(config["q"]) : Array{Int64}([config["q"]])

    seed = config["seed"]
    R = config["R"]
    nev = config["nev"]
    V = config["V"]
    E_window_width = config["E_window_width"]
    num_blas = haskey(config, "num_blas") ? config["num_blas"] : 1
    energy_path = config["energy_path"]
    return Params(l=l, L=L, θ=θ, W=W, E=E, q=q, seed=seed, R=R, nev=nev, V=V,
        E_window_width=E_window_width, num_blas=num_blas, energy_path=energy_path)
end

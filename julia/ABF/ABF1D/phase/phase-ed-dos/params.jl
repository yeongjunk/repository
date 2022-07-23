using Parameters

@with_kw struct Params
    L::Int64 = 100
    θ::Vector{Float64} = [0.25]
    W::Vector{Float64} = [0.]
    E_edges::Vector{Float64} = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0] # Normalized energy, [0, 1]]
    seed::Int64 = 1234
    R::Int64 = 100
    V::Float64 = 1.
    num_blas::Int64 = 1 
    energy_path::String
end

function readconfig(config::Dict)
    L = config["L"]
    θ = Vector(config["th"]) 
    W = Vector(config["W"])
    E_edges = Vector(config["E_edges"])
    seed = config["seed"]
    R = config["R"]
    V = config["V"]
    num_blas = haskey(config, "num_blas") ? config["num_blas"] : 1
    energy_path = config["energy_path"]
    return Params(L=L, θ=θ, W=W, E_edges=E_edges, seed=seed, R=R, V=V,
          num_blas=num_blas, energy_path=energy_path)
end

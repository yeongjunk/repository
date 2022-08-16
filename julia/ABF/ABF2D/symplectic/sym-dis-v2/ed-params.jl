struct Params
    l::Array{Int64}
    L::Int64
    θ::Vector{Float64}
    W::Vector{Float64}
    E::Vector{Float64}
    q::Vector{Int64}

    seed::Int64
    R::Int64
    nev::Int64
    V1::Float64
    V2::Float64
    bw_auto::Bool

    E_bin_width::Int64
    num_blas::Int64
end

function readconfig(config::Dict)
    θ = range(config["th"][1],config["th"][2], length = config["th"][3])
    seed = config["seed"]
    q = isa(config["q"], AbstractArray) ? Array{Int64}(config["q"]) : Array{Int64}([config["q"]])
    W = range(config["W"][1], config["W"][2], length = config["W"][3])
    E = range(config["E"][1], config["E"][2], length = config["E"][3])

    if typeof(config["l"]) <: Integer
        l = [config["l"]];
    else
        l = config["l"]
    end
    num_blas = haskey(config, "num_blas") ? config["num_blas"] : 1

    return Params(l, config["L"], θ, W, E, q, seed, config["R"], config["nev"],
        config["V1"], config["V2"],
        config["bw_auto"], config["E_bin_width"], num_blas)
end

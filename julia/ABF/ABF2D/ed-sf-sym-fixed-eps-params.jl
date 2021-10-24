struct Params
    l::Int64
    L::Int64
    θ::Vector{Float64}
    seed::Int64
    R::Int64
    q::Vector{Int64}
    E_c::Vector{Float64}
    E_del::Float64
    W::Vector{Float64}
    V1::Float64
    V2::Float64
end

function readconfig(config::Dict)
    θ = range(config["th"][1],config["th"][2], length = config["th"][3])
    seed = config["seed"]
    q = isa(config["q"], AbstractArray) ? Array{Int64}(config["q"]) : Array{Int64}([config["q"]])
    W = range(config["W"][1],config["W"][2], length = config["W"][3])
    E_c = range(config["E_c"][1],config["E_c"][2], length = config["E_c"][3])
    return Params(config["l"], config["L"], θ, seed, config["R"], q, E_c, config["E_del"], W, config["V1"], config["V2"])
end

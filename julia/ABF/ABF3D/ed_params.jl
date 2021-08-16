using Parameters

@with_kw struct Params
    set_E_range::Bool = false
    compute_ipr::Bool = true
    L::Int64 = 10
    θ::Vector{Float64} = [0.25]
    seed::Int64
    R::Int64
    q::Vector{Float64}
    E_min::Float64
    E_max::Float64
end

function readconfig(config::Dict)
    if
        set_E_range = Bool(config["set_E_range"])

    compute_ipr = Bool(config["compute_ipr"])

    L = Int64(config["L"])
    θ =  Array{Float64}(range(config["th"][1],config["th"][2], length = config["th"][3]))
    seed = Int64(config["seed"])
    R = Int64(config["R"])


    if set_E_range
        E_min = Float64(config["E_range"][1])
        E_max = Float64(config["E_range"][2])
    else
        E_min = NaN
        E_max = NaN
    end

    if compute_ipr && haskey(config, "q")
        q = Array{Float64}(range(config["q"][1], config["q"][2], length = config["q"][3]))
    elseif compute_ipr && !haskey(config["q"])
        q = Array{Float64}[2.0] # If q is not provided, this is the default
    else
        q = Array{Float64}[]
    end

    return Params(set_E_range, compute_ipr, L, θ, seed, R, q, E_min, E_max)
end

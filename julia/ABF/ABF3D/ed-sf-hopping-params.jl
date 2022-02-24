struct Params
    set_E_range::Bool
    compute_ipr::Bool
    L::Vector{Int64}
    θ::Vector{Float64}
    seed::Int64
    R::Vector{Int64}
    q::Vector{Float64}
    E_min::Float64
    E_max::Float64
    V::Vector{Float64}
end

function readconfig(config::Dict)
    set_E_range = Bool(config["set_E_range"])
    compute_ipr = Bool(config["compute_ipr"])

    L = isa(config["L"], AbstractArray) ? Array{Int64}(config["L"]) : Array{Int64}([config["L"]]);
    θ = range(config["th"][1],config["th"][2], length = config["th"][3])
    seed = config["seed"]
    R = isa(config["R"], AbstractArray) ? Array{Int64}(config["R"]) : Array{Int64}([config["R"]])
    V = isa(config["V"], AbstractArray) ? Array{Float64}(config["V"]) : Array{Float64}([config["V"]])

    if length(R) != length(L)
        error("Array size of R and L should be the same")
    end

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

    return Params(set_E_range, compute_ipr, L, θ, seed, R, q, E_min, E_max, V)
end

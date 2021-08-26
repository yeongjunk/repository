struct Params
    l::Vector{Int64}
    L::Vector{Int64}
    θ::Vector{Float64}
    seed::Int64
    R::Vector{Int64}
    q::Vector{Int64}
    E_min::Float64
    E_max::Float64
end

function readconfig(config::Dict)
    L = isa(config["L"], AbstractArray) ? Array{Int64}(config["L"]) : Array{Int64}([config["L"]])
    l = isa(config["l"], AbstractArray) ? Array{Int64}(config["l"]) : Array{Int64}([config["l"]])
    θ =  range(config["th"][1],config["th"][2], length = config["th"][3])
    seed = config["seed"]
    R = isa(config["R"], AbstractArray) ? Array{Int64}(config["R"]) : Array{Int64}([config["R"]])
    q = isa(config["q"], AbstractArray) ? Array{Int64}(config["q"]) : Array{Int64}([config["q"]])

    if length(R) != length(L)
        error("Array size of R and L should be the same")
    end
    if length(R) != length(L)
        error("Array size of l and L should be the same")
    end

    return Params(l, L, θ, seed, R, q, E_min, E_max)
end

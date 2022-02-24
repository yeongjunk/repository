struct ABF_Params
    L::Int64
    θ::Vector{Float64}
    seed::Int64
    R::Int64
    q::Vector{Float64}

end

function readconfig(config::Dict)
    L = Int64(config["L"])
    θ =  Array{Float64}(range(config["th"][1],config["th"][2], length = config["th"][3]))
    seed = Int64(config["seed"])
    R = Int64(config["R"])
    q = Array{Float64}(range(config["q"][1], config["q"][2], length = config["q"][3]))
    return Params(L, θ, seed, R, q)
end

struct Params
    L::Int64
    θ::Vector{Float64}
    seed::Int64
    R::Int64
    V::Float64
    num_blas::Int64
    E_window::Float64
end

function readconfig(config::Dict)
    θ = range(config["th"][1],config["th"][2], length = config["th"][3])
    seed = config["seed"]
    num_blas = haskey(config, "num_blas") ? config["num_blas"] : 1

    return Params(config["L"], θ, seed, config["R"], 
        config["V"], num_blas, config["E_window"])

end

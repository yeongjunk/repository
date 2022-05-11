struct Params
    L::Int64
    θ::Vector{BigFloat}
    seed::Int64
    R::Int64
    V::BigFloat
    num_blas::Int64
    E_window::BigFloat
    periodic::Bool
    E_c::BigFloat
    nev::Int64
end

function readconfig(config::Dict)
    θ = big.(range(config["th"][1],config["th"][2], length = config["th"][3]))
    seed = config["seed"]
    num_blas = haskey(config, "num_blas") ? config["num_blas"] : 1

    return Params(config["L"], θ, seed, config["R"], 
        big(config["V"]), num_blas, big(config["E_window"]), config["periodic"], big(config["E_c"]), config["nev"])

end

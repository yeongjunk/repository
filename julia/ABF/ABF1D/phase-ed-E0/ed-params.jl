struct Params
    l::Array{Int64}
    L::Int64
    θ::Vector{Float64}
    W::Vector{Float64}
    E::Vector{Float64}
    q::Vector{Float64}

    seed::Int64
    R::Int64
    nev::Int64

    th_ind_range::Vector{Int64}
    W_ind_range::Vector{Int64}
    E_ind_range::Vector{Int64}

    V::Float64
    bw_auto::Bool

    E_bin_width::Int64
    num_blas::Int64
end

function readconfig(config::Dict)
    θ = range(config["th"][1],config["th"][2], length = config["th"][3])
    seed = config["seed"]
    q = isa(config["q"], AbstractArray) ? Array{Float64}(config["q"]) : Array{Float64}([config["q"]])
    W = range(config["W"][1], config["W"][2], length = config["W"][3])
    E = range(config["E"][1], config["E"][2], length = config["E"][3])

    #-----------Energy scan parameters----------#
    if haskey(config, "E_ind_range")
        if length(config["E_ind_range"]) != 2
            error("Energy index range should be a length 2 vector. [start, end]")
        else
            E_ind_range = config["E_ind_range"]
        end
    else
        E_ind_range = [1, config["E_num"]]
    end

    #-----------Theta scan parameters----------#
    if haskey(config, "th_ind_range")
        if length(config["th_ind_range"]) != 2
            error("Theta index range should be a length 2 vector. [start, end]")
        else
            th_ind_range = config["th_ind_range"]
        end
    else
        th_ind_range = [1, length(θ)]
    end

    #-----------W scan parameters----------#
    if haskey(config, "W_ind_range")
        if length(config["W_ind_range"]) != 2
            error("W index range should be a length 2 vector. [start, end]")
        else
            W_ind_range = config["W_ind_range"]
        end
    else
        W_ind_raange = [1, length(W)]
    end

    if typeof(config["l"]) <: Integer
        l = [config["l"]];
    else
        l = config["l"]
    end
    num_blas = haskey(config, "num_blas") ? config["num_blas"] : 1

    return Params(l, config["L"], θ, W, E, q, seed, config["R"], config["nev"],
        th_ind_range, W_ind_range, E_ind_range, config["V"],
        config["bw_auto"], config["E_bin_width"], num_blas)
end

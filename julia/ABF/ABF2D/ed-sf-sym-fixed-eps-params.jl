struct Params
    l::Int64
    L::Int64
    θ::Vector{Float64}
    seed::Int64
    R::Int64
    q::Vector{Int64}
    E_num::Int64
    start_E_ind::Int64
    end_E_ind::Int64
    W::Vector{Float64}
    V1::Float64
    V2::Float64
    nev::Int64
    bw_auto::Bool
    E_min::Float64
    E_max::Float64
    E_bin_width::Int64
end

function readconfig(config::Dict)
    θ = range(config["th"][1],config["th"][2], length = config["th"][3])
    seed = config["seed"]
    q = isa(config["q"], AbstractArray) ? Array{Int64}(config["q"]) : Array{Int64}([config["q"]])
    W = range(config["W"][1], config["W"][2], length = config["W"][3])
    if config["start_E_ind"] < 1 && config["end_E_ind"] > length(E_c)
        error("Invalid Parameter value: invalid energy index")
    end

    #-----Bandwidth parameters-------#
    if haskey(config, "bw_auto")
        bw_auto = config["bw_auto"]
    else
        bw_auto = true
    end

    if bw_auto == false
        E_min, E_max, E_bin_width = Float64(config["E_min"]), Float64(config["E_max"]), Int64(config["E_bin_width"])
    else
        E_min, E_max, E_bin_width = 0., 0., 0
    end
    #-----------Energy scan parameters handling----------#
    if haskey(config, "start_E_ind") && !haskey(config, "end_E_ind")
        start_E_ind = config["start_E_ind"]
        end_E_ind = config["E_num"]
    elseif !haskey(config, "start_E_ind") && haskey(config, "end_E_ind")
        start_E_ind = 1
        end_E_ind = config["end_E_ind"]
    elseif !haskey(config, "start_E_ind") && !haskey(config, "end_E_ind")
        start_E_ind = 1
        end_E_ind = config["E_num"]
    else
        start_E_ind = config["start_E_ind"]
        end_E_ind = config["end_E_ind"]
    end

    return Params(config["l"], config["L"], θ, seed, config["R"], q,
        config["E_num"],
        start_E_ind, end_E_ind, W,
        config["V1"], config["V2"], config["nev"], config["bw_auto"],
        E_min, E_max, E_bin_width)
end

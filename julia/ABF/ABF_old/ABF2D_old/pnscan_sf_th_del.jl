const modulepath =  "/Users/pcs/codes/chain/ABF2D/module"

using CSV, DataFrames
using ArgParse
using JSON
using Random
using LinearAlgebra

# Load custom modules
push!(LOAD_PATH, modulepath)
using ABF2D
using Lattice
using PN

function abf_pn_sf!(l_pj::Lattice2D, l_fe::Lattice2D,  U, W, V, rng)
    H_fe = copy(l_fe.H)

    D = onsite_dis(l_fe, W, rng)
    T = phase_dis(l_fe, V, rng)

    H_fe .+= (D .+ T)
    H_fe .= U'*H_fe*U   #detangle
    l_pj.H .= project(H_fe)

    e, pn = eig_pn!(Hermitian(Array(l_pj.H)))
    df = DataFrame()
    df.E = e; df.PN = pn

    return df
end

function main(ARGS)
    opts = ArgParseSettings(description="Scan and compute pn for all parameters of nu=2 ABF")
    @add_arg_table! opts begin
    "c"
        help = "configuration"
        arg_type = AbstractString
    end
    # Parse the arguments
    popts   = parse_args(opts)
    config  = JSON.parsefile(popts["c"])

    R = config["R"]
    L = config["L"]
    θ = range(config["th_min"], config["th_max"], length = config["th_num"])
    V = range(config["V_min"], config["V_max"], length = config["V_num"])
    V = 10 .^V
    W = config["W"]

    seed = config["seed"]
    rng = MersenneTwister(seed)

    l_fd = Lattice2D(L, L, 2)
    l_fe = Lattice2D(L, L, 2)
    l_pj = Lattice2D(L, L, 1)
    ham_fd!(l_fd, -1, 1)

    println("Scan started")
    @time for (i, θi) in enumerate(θ)
        # Entangle
        U = U_fe(l_fd, θi)
        l_fe.H .= U*l_fd.H*U'
        @time for (j, Vj) in enumerate(V)
            # Filename
            θi_str = replace.(string(round(θi, digits = 3)), '.'=>"-")
            Vj_str = replace.(string(round(Vj, digits = 3)), '.'=>"-")
            fn = "L$(L)_Th$(θi_str)_V$(Vj_str).csv"

            # Realizations
            df = DataFrame(E = Float64[], PN = Float64[], R = Int64[])
            for r in 1:R
                df_r = abf_pn_sf!(l_pj, l_fe, U, W, Vj, rng)
                df_r.R = fill(r, size(df_r, 1))
                append!(df, df_r)
            end
            CSV.write(fn, df)
            println("Saved:"*fn)
        end
    end
    println("Scan started")
end

main(ARGS)

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

function abf_pn!(l_fe::Lattice2D, l_prj::Lattice2D, U,V, rng)
    T = phase_dis(l_fe,V, rng)
    H_fe = copy(l_fe.H)
    H_fe .+= T
    H_fe .= U'*H_fe*U   #detangle
    l_prj.H .= project(H_fe)
    E, PN = eig_pn!(Hermitian(Array(l_prj.H)))
    df = DataFrame()
    df.E = E; df.PN = PN
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
    L = config["L"];
    θ = config["th"];
    W = config["W"];
    V = config["V"];
    seed = config["seed"]

    θ_str = replace.(string(θ), '.'=>"")
    V = config["V"];

    fn = "L$(L)_Th$(θ_str)"
    rng = MersenneTwister(seed)

    l = Lattice2D(L, L, 2)
    ham_fd!(l, -1, 1)
    U = U_fe(l, θ)
    l.H .= U*l.H*U'
    l_prj = Lattice2D(L,L,1)
    println("start")
    @time for r in 1:R
        fn_r = fn*"_R$(r).csv"
        df = abf_pn!(l, l_prj,U, V, rng)
        CSV.write(fn_r, df)
    end
end

main(ARGS)

# Compute the eigenvalue using Arnoldi method & shift-and-inverse technique.
using ArgParse, JSON
using LinearAlgebra, Arpack, LinearMaps
using SparseArrays
using Random
using DelimitedFiles

using ABF
using Lattice

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

    L = Int64(config["L"])
    θ = Float64(config["th"])
    seed = Int64(config["seed"])
    E = Float64(config["E"])
    nev = Int64(config["nev"])
    R = Int64(config["R"])
    W = 1.;

    fn = "L_$(L)_Th_$(θ)_E_$(E)"
    rng = MersenneTwister(seed)

    l = Lattice3D(L, L, L, 2)
    l_p = Lattice3D(L, L, L, 1)

    H = ham_fd(l, -1, 1)
    U = U_fe(l, θ)
    H .= U*H*U'
    data = Array{Float64}(undef, L^3+1, nev)
    for r in 1:R
        D = spdiagm(0 => W*(rand(rng, size(H,1)) .- 0.5))

        H_p = project(dropzeros(U'*(H + D)*U))
        H_p .= sparse(Symmetric(H_p))
        @time data[1,:], data[2:end,:] = eigs(H_p, nev = nev, sigma = E)
        open(fn*"_r_$(r).csv", "w") do io
            writedlm(io, data, ',')
        end
    end
end

main(ARGS)

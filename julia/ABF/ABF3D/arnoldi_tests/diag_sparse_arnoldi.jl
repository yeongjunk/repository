# Compute the eigenvalue using Arnoldi method & shift-and-inverse technique.

using ArgParse, JSON
using LinearAlgebra, ArnoldiMethod, LinearMaps
using SparseArrays
using Random
using JLD2

using ABF
using Lattice

# The linear map for the inverse matrix
function construct_linear_map(A)
    F = lu(A)
    LinearMap{eltype(A)}((y, x) -> ldiv!(y, F, x), size(A,1), ismutating=true)
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
        D = spdiagm(0 => W*(rand(size(H,1)) .- 0.5))

        H_p = project(dropzeros(U'*(H + D)*U))
        H_p .= sparse(Symmetric(H_p))
        droptol!(H_p, 1E-12)
        σI = spdiagm(0 => fill(E, L^3))
        @time decomp, = partialschur(construct_linear_map(H_p-σI), nev=nev, tol=1e-6, restarts=100, which=LM())
        Es_shift_inv, ψs = partialeigen(decomp)
        # Eigenvalues have to be inverted to find the smallest eigenvalues of the non-inverted problem.
        Es = (1 ./ Es_shift_inv) .+ E

        data[1,:] = Es
        data[2:end,:] = ψs
        JLD2.save(fn*"_r_$(r).jld2", Dict("data" => data))
    end
end

main(ARGS)

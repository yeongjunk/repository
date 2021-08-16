using LinearAlgebra, SparseArrays
using Random
using DelimitedFiles
using ArgParse, JSON

# Custom modules
using ABF
using Lattice
using PN

include("./ed_params.jl") # read parameters from configuration file
function abf2d_scan(p::Params)
    rng = MersenneTwister(p.seed)
    l = Lattice2D(p.L, p.L, 2)
    l_prj = Lattice2D(p.L, p.L, 1)
    for i in 1:length(p.θ)
        fn = "L_$(p.L)_Th_$(i)_R_$(p.R)" #File name

        H, U = ham_fe(l, -1, 1, p.θ[i]) # Fully entangled hamiltonian

        @time for r in 1:p.R
            # Add disorder & detangle & project
            D = spdiagm(0 => rand(rng, size(H,1)) .- 0.5)
            H_prj = project(U'*(H + D)*U)
            droptol!(H_prj, 1E-12)

            if p.set_E_range
                e, psi = eigen!(Symmetric(Matrix(H_prj)), p.E_min, p.E_max)
            else
                e, psi = eigen!(Symmetric(Matrix(H_prj)))
            end
            e = transpose(e)
            data = vcat(e, psi)

            open(fn*"_r_$(r).csv", "w") do io
                writedlm(io, data, ',')
            end

        end
    end
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
    p = readconfig(config)
    abf2d_scan(p)
end

main(ARGS)

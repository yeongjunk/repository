using LinearAlgebra, SparseArrays
using Random
using DataFrames, CSV
using ArgParse, JSON

# Custom modules
using ABF
using Lattice
using PN

include("./ed_params.jl") # read parameters from configuration file
function abf2d_scan(p::Params)
    q_str = string.(p.q)
    rng = MersenneTwister(p.seed)
    l = Lattice2D(p.L, p.L, 2)
    l_prj = Lattice3D(p.L, p.L, p.L, 1)
    for i in 1:length(p.θ)
        fn = "L$(p.L)_Th$(i)_R$(p.R).csv" #File name

        H, U = ham_fe(l, -1, 1, p.θ[i]) # Fully entangled hamiltonian
        df = DataFrame(E = Float64[], r = Int64[])
        if p.compute_ipr
            for j in 1:length(p.q)
                insertcols!(df, q_str[j] => Float64[])
            end
        end
        @time for r in 1:p.R
            # Add disorder & detangle & project
            D = spdiagm(0 => rand(rng, size(H,1)) .- 0.5)
            H_prj = project(U'*(H + D)*U)
            droptol!(H_prj, 1E-12)

            if p.compute_ipr && p.set_E_range
                e, psi = eigen!(Symmetric(Matrix(H_prj)), p.E_min, p.E_max)
            elseif p.compute_ipr && !p.set_E_range
                e, psi = eigen!(Symmetric(Matrix(H_prj)))

            elseif !p.compute_ipr && p.set_E_range
                e = eigvals!(Symmetric(Matrix(H_prj)), p.E_min, p.E_max)
            elseif !p.compute_ipr && !p.set_E_range
                e = eigvals!(Symmetric(Matrix(H_prj)))
            end

            df_temp = DataFrame(E = e, r = fill(r, length(e)))
            if p.compute_ipr
                for j in 1:length(p.q)
                    insertcols!(df_temp, q_str[j] => round.(compute_iprs(psi, q = p.q[j]), sigdigits = 6))
                end
            end

            append!(df, df_temp)
        end
        CSV.write(fn, df)
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

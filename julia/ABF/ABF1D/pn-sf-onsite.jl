using LinearAlgebra, SparseArrays
using Random
using DataFrames, CSV
using ArgParse, JSON

# Custom modules
using ABF
using Lattice
using PN

include("./pn-sf-onsite-params.jl") # read parameters from configuration file
function abf1d_scan(p::Params)
    q_str = "q=".*string.(p.q)
    rng = MersenneTwister(p.seed)
    for (i, Li) in enumerate(p.L)
        ltc = Lattice1D(Li, 2)
        for (j, θj) in enumerate(p.θ)
            fn = "d1abf-sf-onsite-pn-L$(Li)-th$(round(θj, digits = 3))-R$(p.R[i]).csv" #File name
            H, U = ham_fe(ltc, -1., 1., θj) # Fully entangled hamiltonian
            df = DataFrame(E = Float64[], r = Int64[])
            if p.compute_ipr
                for j in 1:length(p.q)
                    insertcols!(df, q_str[j] => Float64[])
                end
            end
            @time for r in 1:p.R[i]
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

                df_temp = DataFrame(E = round.(e, sigdigits = 10), r = fill(r, length(e)))
                if p.compute_ipr
                    for j in 1:length(p.q)
                        insertcols!(df_temp, q_str[j] => round.(compute_iprs(psi, q = p.q[j]), sigdigits = 10))
                    end
                end
                append!(df, df_temp)
            end
            CSV.write(fn, df)
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
    abf1d_scan(p)
end

main(ARGS)

using LinearAlgebra, SparseArrays
using Random
using DataFrames, CSV
using ArgParse, JSON

# Custom modules
using ABF
using Lattice
using PN

function phase_dis(H; V = 1., rng = nothing)
    D = convert.(ComplexF64, H)
    if rng == nothing
        rng = MersenneTwister()
    end
    rows = rowvals(D)
    vals = nonzeros(D)
    m, n = size(D)
    for j = 1:n
       for i in nzrange(D, j)
          row = rows[i]
          # println("$row ,", "$j")
          if row > j
              vals[i] = im*V*vals[i]*(rand(rng) .- 0.5)
          elseif row <= j
              vals[i] = 0.
          end
       end
    end
    return D + D'
end

function hop_dis(H; V = 1., rng = nothing)
    D = convert.(ComplexF64, H)
    if rng == nothing
        rng = MersenneTwister()
    end
    rows = rowvals(D)
    vals = nonzeros(D)
    m, n = size(D)
    for j = 1:n
       for i in nzrange(D, j)
          row = rows[i]
          # println("$row ,", "$j")
          if row > j
              vals[i] = V*vals[i]*(rand(rng) .- 0.5)
          elseif row <= j
              vals[i] = 0.
          end
       end
    end
    return D + D'
end


include("./pn-sf-phase-params.jl") # read parameters from configuration file
function abf1d_scan(p::Params)
    q_str = "q=".*string.(p.q)
    rng = MersenneTwister(p.seed)
    for (i, Li) in enumerate(p.L)
        ltc = Lattice1D(Li, 2)
        for (j, θj) in enumerate(p.θ)
            for (k, Vk) in enumerate(p.V)
                fn = "d1abf-sf-phase-pn-L$(Li)-th$(round(θj, digits = 3))-V$(Vk)-R$(p.R[i]).csv" #File name
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
                    D_phase = phase_dis(H, V = Vk)
                    H_prj = project(U'*(H + D + D_phase)*U)
                    droptol!(H_prj, 1E-12)

                    if p.compute_ipr && p.set_E_range
                        e, psi = eigen!(Hermitian(Matrix(H_prj)), p.E_min, p.E_max)
                    elseif p.compute_ipr && !p.set_E_range
                        e, psi = eigen!(Hermitian(Matrix(H_prj)))

                    elseif !p.compute_ipr && p.set_E_range
                        e = eigvals!(Hermitian(Matrix(H_prj)), p.E_min, p.E_max)
                    elseif !p.compute_ipr && !p.set_E_range
                        e = eigvals!(Hermitian(Matrix(H_prj)))
                    end

                    df_temp = DataFrame(E = e, r = fill(r, length(e)))
                    if p.compute_ipr
                        for j in 1:length(p.q)
                            insertcols!(df_temp, q_str[j] => round.(compute_iprs(psi, q = p.q[j]), digits = 10))
                        end
                    end
                    append!(df, df_temp)
                end
                CSV.write(fn, df)
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
    abf1d_scan(p)
end

main(ARGS)

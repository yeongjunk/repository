using LinearAlgebra, SparseArrays, Arpack
using Random
using DataFrames, CSV
using ArgParse, JSON

# Custom modules
using ABF
using Lattice
using PN

function box_inds(ltc, b)
    L = ltc.L
    box_pts = b^3
    box_num = L^3 ÷ b^3

    inds = Array{Int64}(undef, box_pts, box_num)
    box_count = 1
    for x in 1:b:(L-b+1), y in 1:b:(L-b+1), z in 1:b:(L-b+1)
        for l in 0:b-1, m in 0:b-1, n in 0:b-1
            inds[b^2*l + b*m + n + 1, box_count] = index(ltc, (x + l, y + m, z + n, 1))
        end
        box_count += 1
    end
    return inds
end

function compute_box_iprs(ltc, psi, idx::Array{Int64}; q = 2)
    p = Array{Float64}(undef, size(idx, 2), size(psi, 2))
    for i in 1:size(psi, 2)
        for j in 1:size(idx, 2)
            @views p[j, i] = sum(abs2,  psi[idx[:, j], i])
        end
    end
    return compute_iprs2(p, q = q)
end


function compute_box_iprs(ltc, psi, b::Int64; q = 2)
    idx = box_inds(ltc, b)
    p = Array{Float64}(undef, size(idx, 2), size(psi, 2))
    for i in 1:size(psi, 2)
        for j in 1:size(idx, 2)
            @views p[j, i] = sum(abs2, psi[idx[:, j], i])
        end
    end
    return compute_iprs2(p)
end


include("./ed-sf-onsite-fixed-eps-params.jl") # read parameters from configuration file
function abf3d_scan(p::Params)
    q_str = "q".*string.(p.q)
    rng = MersenneTwister(p.seed)
    for i in 1:length(p.L)
        ltc = Lattice3D(p.L[i], p.L[i], p.L[i], 2)
        ltc_p = Lattice3D(p.L[i], p.L[i], p.L[i], 1)
        boxidx = box_inds(ltc_p, p.l[i])
        for j in 1:length(p.θ)
            fn = "L$(p.L[i])_Th$(i).csv" #File name

            H, U = ham_fe(ltc, -2, 0, p.θ[j]) # Fully entangled hamiltonian
            df = DataFrame(E = Float64[], r = Int64[])
            for k in 1:length(p.q)
                insertcols!(df, q_str[k] => Float64[])
            end
            r = 1
            @time while size(df, 1) <= p.R[i]
                # Add disorder & detangle & project
                D = Diagonal(rand(rng, size(H,1)) .- 0.5)
                H_prj = project(U'*(H + D)*U)
                droptol!(H_prj, 1E-12)

                e, psi = eigs(Symmetric(H_prj), nev = p.L[i]^3÷500, sigma = 0.)
                idx = findall(x -> E_min < x && x < E_max, e)
                @views df_temp = DataFrame(E = round.(e[idx], sigdigits = 12), r = fill(r, length(e)))
                for k in 1:length(p.q)
                    @views insertcols!(df_temp, q_str[k] => round.(compute_box_iprs(ltc_p, psi[:, idx], boxidx, q = p.q[k]), sigdigits = 12))
                end
                append!(df, df_temp)
                r += 1
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
    abf3d_scan(p)
end

main(ARGS)

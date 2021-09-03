using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using DataFrames, CSV
using ArgParse, JSON
# Custom modules
using ABF
using Lattice
using PN
function construct_linear_map(A)
    F = factorize(A)
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), size(A, 1), ismutating = true)
end

function ldiv2!(y, F, x)
    y .= F\x
end

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
          #println("$row ,", "$j")
          if row > j
              vals[i] = im*V*vals[i]*(rand(rng) .- 0.5)
          elseif row <= j
              vals[i] = 0.
          end
       end
    end
    return D + D'
end

function box_inds(ltc, b)
    L = ltc.M
    box_pts = b^2
    box_num = L^2 ÷ b^2

    inds = Array{Int64}(undef, box_pts, box_num)
    box_count = 1
    for y in 1:b:(L-b+1), z in 1:b:(L-b+1)
        for m in 0:b-1, n in 0:b-1
            inds[b*m + n + 1, box_count] = index(ltc, (y + m, z + n, 1))
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
    return compute_iprs2(p, q = q)
end


include("./ed-sf-onsite-fixed-eps-params.jl") # read parameters from configuration file
function abf3d_scan(p::Params)
    q_str = "q".*string.(p.q)
    rng = MersenneTwister(p.seed)
    ltc = Lattice2D(p.L, p.L, 2)
    ltc_p = Lattice2D(p.L, p.L, 1)
    boxidx = box_inds(ltc_p, p.l)
    for j in 1:length(p.θ)
        fn = "L$(p.L)_Th$(j)" #File name
        H, U = ham_fe(ltc, -2, 0, p.θ[j]) # Fully entangled hamiltonian
        df = DataFrame(E = Float64[], r = Int64[])
        for k in 1:length(p.q)
            insertcols!(df, q_str[k] => Float64[])
        end

        for jj in 1:length(p.W), jjj in 1:length(p.E_c)
            r = 1
            @time while size(df, 1) <= p.R
                # Add disorder & detangle & project
                D = phase_dis(H, rng = rng) .+ p.W[jj]*Diagonal(rand(rng, size(H,1)) .- 0.5)
                @views H_prj = project(U'*D*U)
                droptol!(H_prj, 1E-12)
                e_inv, psi, info = eigsolve(construct_linear_map(Hermitian(H_prj .- p.E_c[jjj]*I(size(H_prj, 1)))), size(H_prj, 1), div(p.L^2, 400), :LM, ishermitian = true, krylovdim = max(30, 2*div(p.L^2, 400)+1));
                e = 1 ./ real.(e_inv) .+ p.E_c[jjj]
                psi = reduce(hcat, psi)
                idx = findall(x -> (p.E_c[jjj] - p.E_del) < x && x < (p.E_c[jjj] + p.E_del), e)
                @views df_temp = DataFrame(E = round.(e[idx], sigdigits = 12), r = fill(r, length(idx)))
                for k in 1:length(p.q)
                    @views insertcols!(df_temp, q_str[k] => round.(compute_box_iprs(ltc_p, psi[:, idx], boxidx, q = p.q[k]), sigdigits = 12))
                end
                append!(df, df_temp)
                r += 1
            end
            CSV.write(fn*"_W$(jj)_E$(jjj).csv", df)
            df = DataFrame(E = Float64[], r = Int64[])
            for k in 1:length(p.q)
                insertcols!(df, q_str[k] => Float64[])
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
    abf3d_scan(p)
end

main(ARGS)

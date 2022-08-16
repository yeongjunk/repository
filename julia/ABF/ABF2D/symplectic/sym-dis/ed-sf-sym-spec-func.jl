using LinearAlgebra, SparseArrays, Arpack
using LinearMaps
using Random
using DataFrames, CSV
using ArgParse, JSON
# Custom modules
using ABFSym
using Lattice
using PN
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
include("./ed-sf-sym-fixed-eps-params.jl") # read parameters from configuration file

function construct_linear_map(A)
    F = factorize(A)
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), size(A, 1), ismutating = true)
end

function ldiv2!(y, F, x)
    y .= F\x
end

function symplectic_coupling!(ltc, t, I, J, val, site1, site2, u, v, V)
    push!(I, index(ltc, (site1[1], site1[2], u)))
    push!(J, index(ltc, (site2[1], site2[2], v)))
    push!(val, t*V[1, 1])

    push!(I, index(ltc, (site1[1], site1[2], u)))
    push!(J, index(ltc, (site2[1], site2[2], v + 2)))
    push!(val, t*V[1, 2])

    push!(I, index(ltc, (site1[1], site1[2], u + 2)))
    push!(J, index(ltc, (site2[1], site2[2], v)))
    push!(val, t*V[2, 1])

    push!(I, index(ltc, (site1[1], site1[2], u + 2)))
    push!(J, index(ltc, (site2[1], site2[2], v + 2)))
    push!(val,  t*V[2 ,2])
end


function symp_dis(V1, V2; rng = Random.GLOBAL_RNG)
    V2_dis = V2*(rand(rng) .- 0.5)
    ϕ_dis = 2pi*rand(rng)
    return [V1 exp(-im*ϕ_dis)*V2_dis; -exp(im*ϕ_dis)*V2_dis V1]
end

function makesym2d(ltc, H, V1, V2; rng = Random.GLOBAL_RNG)
    t0 = Float64[]
    t1 = Float64[]
    t2 = Float64[]
    t3 = Float64[]
    t4 = Float64[]

    n, m = 1, 1
    for i in 1:2, j in 1:2
        push!(t0, H[index(ltc, (m, n, i)), index(ltc, (m, n, j))] )
        push!(t1, H[index(ltc, (m, n, i)), index(ltc, (m, n + 1, j))])
        push!(t2, H[index(ltc, (m, n, i)), index(ltc, (m + 1, n, j))])
        push!(t3, H[index(ltc, (m, n, i)), index(ltc, (m + 1, n + 1, j))])
        push!(t4, H[index(ltc, (m + 1, n, i)), index(ltc, (m, n + 1, j))])
    end
    I, J, val = findnz(sparse(Diagonal(H)))
    val = Vector(val)
    for m in 1:ltc.M, n in 1:ltc.N
        for i in 1:2, j in 1:2
            if i != j
                symplectic_coupling!(ltc, t0[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m, n), i, j, [1 0; 0 1])
            end
            Vx = symp_dis(V1, V2, rng = rng)
            symplectic_coupling!(ltc, t1[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m, n + 1), i, j, Vx)
            symplectic_coupling!(ltc, t1[2(i-1) + (j-1) + 1], I, J, val, (m, n + 1), (m, n), j, i, Vx')
            Vy = symp_dis(V1, V2, rng = rng)
            symplectic_coupling!(ltc, t2[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m + 1, n), i, j, Vy)
            symplectic_coupling!(ltc, t2[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n), (m, n), j, i, Vy')
            Vd1 = symp_dis(V1, V2, rng = rng)
            symplectic_coupling!(ltc, t3[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m + 1, n + 1), i, j, Vd1)
            symplectic_coupling!(ltc, t3[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n + 1), (m, n), j, i, Vd1')
            Vd2 = symp_dis(V1, V2, rng = rng)
            symplectic_coupling!(ltc, t4[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n), (m, n + 1), i, j, Vd2)
            symplectic_coupling!(ltc, t4[2(i-1) + (j-1) + 1], I, J, val, (m, n + 1), (m + 1, n), j, i, Vd2')
        end
    end
    return sparse(I, J, val, size(H, 1), size(H, 2))
end
function abf3d_scan(p::Params)
    rng = MersenneTwister(p.seed)
    ltc = Lattice2D(p.L, p.L, 4)
    ltc_p = Lattice2D(p.L, p.L, 2)
    for j in 1:length(p.θ)
        fn = "spec_L$(p.L)_Th$(j)" #File name
        H, U = ham_fe(ltc, -2, 0, p.θ[j]) # Fully entangled hamiltonian
        H = convert.(ComplexF64, H)
        df = DataFrame(E = Float64[], r = Int64[])
        for jj in 1:length(p.W), jjj in 1:length(p.E_c)
            r = 1
            @time while size(df, 1) <= p.R
                # Add disorder & detangle & project
                H_dis = makesym2d(ltc, H, p.V1, p.V2) + p.W[jj]*Diagonal(rand(rng, size(H,1)) .- 0.5)
                @views H_prj = project(U'*H_dis*U)
                droptol!(H_prj, 1E-12)
                e_inv, _ = eigs(construct_linear_map(500. * Hermitian(H_prj .- p.E_c[jjj]*I(size(H_prj, 1)))), nev = div(p.L^2, 100), which = :LM, ritzvec=false);
                e = 1 ./ (500. * real.(e_inv)).+ p.E_c[jjj]
                idx = findall(x -> (p.E_c[jjj] - p.E_del) < x && x < (p.E_c[jjj] + p.E_del), e)
                @views df_temp = DataFrame(E = round.(e[idx], sigdigits = 12), r = fill(r, length(idx)))
                append!(df, df_temp)
                r += 1
            end
            CSV.write(fn*"_W$(jj)_E$(jjj).csv", df)
            df = DataFrame(E = Float64[], r = Int64[])
        end
    end
end

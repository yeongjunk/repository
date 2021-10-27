
using LinearAlgebra, SparseArrays, ArnoldiMethod
using LinearMaps
using Random
using DataFrames, CSV
using ArgParse, JSON
# Custom modules
using ABFSym
using Lattice
using PN

LinearAlgebra.BLAS.set_num_threads(1)
include("./ed-sf-sym-fixed-eps-params.jl") # read parameters from configuration file

function construct_linear_map(A)
    F = factorize(A)
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), size(A, 1), ismutating = true)
end

function ldiv2!(y, F, x)
    y .= F\x
end

function box_inds(ltc, b)
    L = ltc.M
    U = ltc.U
    box_pts = b^2
    box_num = U*L^2 ÷ b^2

    inds = Array{Int64}(undef, box_pts, box_num)
    box_count = 1
    for y in 1:b:(L-b+1), z in 1:b:(L-b+1), u in 1:U
        for m in 0:b-1, n in 0:b-1
            inds[b*m + n + 1, box_count] = index(ltc, (y + m, z + n, u))
        end
        box_count += 1
    end
    return inds
end

# @doc"""
# ltc: 2D Lattice
# psi: wavefunction vector on the lattice support
# idx: box indices
# q: GIPR exponent
# """
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

function estimate_bw(p, θ, W)
    ltc = Lattice2D(20, 20, 4)
    ltc_p = Lattice2D(20, 20, 2)
    rng = MersenneTwister()
    H, U = ham_fe(ltc, -2, 0, θ) # Fully entangled hamiltonian
    H = convert.(ComplexF64, H)
    H_dis = makesym2d(ltc, H, p.V1, p.V2, rng = rng)
    D = W*Diagonal(rand(rng, size(H,1)) .- 0.5)
    @views H_prj = project(U'*(H_dis + D)*U)
    vals = eigvals!(Hermitian(Matrix(H_prj)))
    return 2*maximum(vals)
end
function abf3d_scan(p::Params)
    q_str = "q".*string.(p.q)
    rng = MersenneTwister(p.seed)
    ltc = Lattice2D(p.L, p.L, 4)
    ltc_p = Lattice2D(p.L, p.L, 2)
    boxidx = box_inds(ltc_p, p.l)
    DF_store = Array{DataFrame}(undef, length(p.θ), length(p.W), length(p.E_c))
    FN_store = Array{String}(undef, length(p.θ), length(p.W), length(p.E_c))
    println("Number of threads: $(Threads.nthreads()). Start scanning")
    @Threads.threads for j in 1:length(p.θ)
        fn = "L$(p.L)_Th$(j)" #File name
        H, U = ham_fe(ltc, -2, 0, p.θ[j]) # Fully entangled hamiltonian
        H = convert.(ComplexF64, H)
        df = DataFrame(E = Float64[], r = Int64[])
        for k in 1:length(p.q)
            insertcols!(df, q_str[k] => Float64[])
        end
        for jj in 1:length(p.W)
            BW = estimate_bw(p, p.θ[j], p.W[jj])
            println("Bandwidth estimated: $(BW)")
            E_c = range(0.0001, BW/2*0.95, length = length(p.E_c))
            E_del = (E_c[2] - E_c[1])/8
            for jjj in 1:length(E_c)
                r = 1
                @time while size(df, 1) <= p.R
                    # Add disorder & detangle & project
                    H_dis = makesym2d(ltc, H, p.V1, p.V2, rng = rng)
                    D = p.W[jj]*Diagonal(rand(rng, size(H,1)) .- 0.5)
                    @views H_prj = project(U'*(H_dis + D)*U)
                    droptol!(H_prj, 1E-12)
                    decomp, = partialschur(construct_linear_map(500. * Hermitian(H_prj .- E_c[jjj]*I(size(H_prj, 1)))), nev = div(p.L^2, 100), tol=1e-6, restarts=100, which = LM())
                    e_inv, psi = partialeigen(decomp)
                    e = 1 ./ (500. * real.(e_inv)) .+ E_c[jjj]
                    idx = findall(x -> (E_c[jjj] - E_del) < x && x < (E_c[jjj] + E_del), e)
                    @views df_temp = DataFrame(E = round.(e[idx], sigdigits = 12), r = fill(r, length(idx)))
                    for k in 1:length(p.q)
                        @views insertcols!(df_temp, q_str[k] => round.(compute_box_iprs(ltc_p, psi[:, idx], boxidx, q = p.q[k]), sigdigits = 12))
                    end
                    append!(df, df_temp)
                    r += 1
                end
                DF_store[j, jj, jjj] = df
                FN_store[j, jj, jjj] = fn*"_W$(jj)_E$(jjj).csv"
                df = DataFrame(E = Float64[], r = Int64[])
                for k in 1:length(p.q)
                    insertcols!(df, q_str[k] => Float64[])
                end
            end
        end
    end
    for i in eachindex(DF_store) #save
        CSV.write(FN_store[i], DF_store[i])
    end
end

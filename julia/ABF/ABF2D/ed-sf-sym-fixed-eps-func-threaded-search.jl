
using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using DataFrames, CSV, DelimitedFiles
using ArgParse, JSON
# Custom modules
using ABFSym
using Lattice
using PN
using Glob
include("./ed-sf-sym-fixed-eps-params.jl") # read parameters from configuration file

function construct_linear_map(A)
    F = factorize(A)
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), size(A, 1), ismutating = true, ishermitian = true)
end

function ldiv2!(y, F, x)
    y .= F\x
end

function dis(N, W, rng)
    arr = Vector{Float64}(undef, 4N)
    d_arr = W*(rand(rng, 2N) .- 0.5)
    for i in 1:N
        arr[4(i-1) + 1:4(i-1)+1 + 1] = d_arr[2(i-1) + 1:2(i-1) + 1 + 1]
        arr[4(i-1) + 1 + 2:4(i-1) + 1 + 3] = d_arr[2(i-1) + 1:2(i-1) + 1 + 1]
    end
    return arr
end

function box_inds(ltc, b)
    L = ltc.M
    U = ltc.U
    box_pts = U*b^2
    box_num = L^2 ÷ b^2

    inds = Array{Int64}(undef, box_pts, box_num)
    box_count = 1
    for y in 1:b:(L-b+1), z in 1:b:(L-b+1)
        for m in 0:b-1, n in 0:b-1, u in 1:U
            inds[2b*m + 2*n + u, box_count] = index(ltc, (y + m, z + n, u))
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
    I = Int64[];
    J = Int64[];
    val = ComplexF64[];
    Id = LinearAlgebra.I(2)
    for m in 1:ltc.M, n in 1:ltc.N
        for i in 1:2, j in 1:2

            Vx = symp_dis(V1, V2, rng = rng)
            Vy = symp_dis(V1, V2, rng = rng)
            Vd1 = symp_dis(V1, V2, rng = rng)
            Vd2 = symp_dis(V1, V2, rng = rng)
            symplectic_coupling!(ltc, t0[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m, n), i, j, Id)
            symplectic_coupling!(ltc, t1[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m, n + 1), i, j, Vx)
            symplectic_coupling!(ltc, t1[2(i-1) + (j-1) + 1], I, J, val, (m, n + 1), (m, n), j, i, Vx')
            symplectic_coupling!(ltc, t2[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m + 1, n), i, j, Vy)
            symplectic_coupling!(ltc, t2[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n), (m, n), j, i, Vy')
            symplectic_coupling!(ltc, t3[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m + 1, n + 1), i, j, Vd1)
            symplectic_coupling!(ltc, t3[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n + 1), (m, n), j, i, Vd1')
            symplectic_coupling!(ltc, t4[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n), (m, n + 1), i, j, Vd2)
            symplectic_coupling!(ltc, t4[2(i-1) + (j-1) + 1], I, J, val, (m, n + 1), (m + 1, n), j, i, Vd2')
        end
    end
    return sparse(I, J, val, size(H, 1), size(H, 2))
end

function estimate_bw(p, θ, W, L, rng)
    ltc = Lattice2D(L, L, 4)
    H, U = ham_fe(ltc, -2, 0, θ) # Fully entangled hamiltonian
    H = convert.(ComplexF64, H)
    H_dis = makesym2d(ltc, H, p.V1, p.V2, rng = rng)
    D = Diagonal(dis(L^2, W, rng))
    @views H_prj = project(U'*(H_dis + D)*U)
    vals, psi, info = eigsolve(Hermitian(H_prj), size(H_prj, 1), 1, :LM, ishermitian = true, krylovdim = 30)
    return 2*abs(maximum(vals))
end

# @doc"""
# p: Params
# θ: angle
# W: disorder strength
# L: system size for estimating bandwidth.
# E_crop: crop the end of the spectrum. Recommend 0.9
# E_bin_width: portion given by 1/E_bin_width of the small energy window around energy center.
# rng: random number generator
# """
function auto_energy_params(p, θ, W, L, E_crop, E_bin_width, rng)
    BW = estimate_bw(p, θ, W, L, rng)
    E_c = range(0.0001, BW/2*E_crop, length = length(p.E))
    E_del = (E_c[2] - E_c[1])/E_bin_width
    return BW, E_c, E_del
end

function energy_param_generator(p, θ, W, L, E_crop, E_bin_width, rng)
    if p.bw_auto
        BW, E_c, E_del = auto_energy_params(p, θ, W, L, E_crop, E_bin_width, rng)
    else
        if length(p.E) != 1
            BW = p.E[end] - p.E[1]
            E_c = p.E
            E_del = (E_c[2] - E_c[1])/p.E_bin_width
        elseif length(p.E) == 1
            BW = estimate_bw(p, θ, W, L, rng)
            E_c = p.E[1] + 1E-12
            E_del = BW/p.E_bin_width
        end
    end
    return BW, E_c, E_del
end

function print_info(nt, BW, E_c, E_del, p)
    println("Number of BLAS threads: ",LinearAlgebra.BLAS.get_num_threads())
    println("Number of threads: $(nt)")
    println("Realization / nev: ", p.R/p.nev)
    println("Auto bandwidth settings: ", p.bw_auto)
    println("Bandwidth estimated: $(BW)")
    println("Energy centers: ", E_c)
    println("Energy bin width: $(E_del)")
    println("Thetas: ", p.θ)
end

function generate_col_names(p::Params)
    return ["l"*string(li)*"_q"*string(qi) for li in p.l, qi in p.q]
end

function missing_idx_finder(p::Params)
    idcs = Tuple{Int64, Int64, Int64}[]
    lists = glob("L*_Th*_W*_E*.csv")
    miss_idx = (0, 0, 0)
    for i in 1:length(lists)
        fn = split(split(lists[i], ".")[1], "_")
        x = parse(Int64, replace(fn[2], "Th"=>""))
        y = parse(Int64, replace(fn[3], "W"=>""))
        z = parse(Int64, replace(fn[4], "E"=>""))
        push!(idcs, (x, y, z))
    end
    for l in 1:length(p.θ), m in 1:length(p.W), n in 1:length(p.E)
        if length(findall(x -> x == (l, m, n), idcs)) == 0
            miss_idx = (l, m ,n)
            break
        end
    end
    return miss_idx
end



function abf3d_scan(p::Params)
    col_str = generate_col_names(p)
    nt = Threads.nthreads()
    rng = [MersenneTwister(p.seed + i) for i in 1:nt]
    ltc = Lattice2D(p.L, p.L, 4)
    ltc_p = Lattice2D(p.L, p.L, 2)
    boxidx = [box_inds(ltc_p, li) for li in p.l]
    tag = rand(UInt64)
    while true
        @label Search
        j, jj, jjj = missing_idx_finder(p)
        if (j, jj, jjj) == (0, 0, 0)
            println("It seems that scan is finished. The scan will be terminated.")
            break
        end
        CSV.write("L$(p.L)_Th$(j)_W$(jj)_E$(jjj)_temp_$(tag).csv", DataFrame())
        scan_overlap = glob("L$(p.L)_Th$(j)_W$(jj)_E$(jjj)_temp_*.csv")
        sleep(3rand())
        if length(scan_overlap) != 1
            println("Overlap detected")
            rm("L$(p.L)_Th$(j)_W$(jj)_E$(jjj)_temp_$(tag).csv")
            @goto Search
        end
        println("start scanning at index: ", (j, jj, jjj))
        H, U = ham_fe(ltc, -2, 0, p.θ[j]) # Fully entangled hamiltonian
        BW, E_c, E_del = energy_param_generator(p, p.θ[j], p.W[jj], 50, 0.9, 8, rng[1])
        # print_info(nt, BW, E_c, E_del, p)
        df = [DataFrame(E = Float64[], r = Int64[]) for i in 1:nt]
        for t in 1:nt, k in 1:length(p.l), kk in 1:length(p.q)
            insertcols!(df[t], col_str[k, kk] => Float64[])
        end
        @Threads.threads for r in 1:p.R÷p.nev# Realizations
            er = true
            er_num = 0
            while er
                er_num == 5 && error("5 attempts faild.")
                try
                    x = Threads.threadid()
                    H_dis = makesym2d(ltc, H, p.V1, p.V2, rng = rng[x])
                    D = Diagonal(dis(p.L^2, p.W[jj], rng[x]))
                    @views H_prj = project(U'*(H_dis + D)*U)
                    droptol!(H_prj, 1E-12)
                    lmap = construct_linear_map(Hermitian(p.L^2*(H_prj - E_c[jjj]*I(size(H_prj, 1)))))
                    e_inv, psi, info = eigsolve(lmap, size(H_prj, 1), p.nev, :LM, ishermitian = true, krylovdim = max(30, 2p.nev + 1));
                    e = 1 ./ (p.L^2*real.(e_inv)) .+ E_c[jjj]
                    psi = reduce(hcat, psi)
                    #Crop energies that are outside the energy bins
                    idx = findall(x -> (E_c[jjj] - E_del) <  x < (E_c[jjj] + E_del), e)
                    @views df_temp = DataFrame(E = round.(e[idx], sigdigits = 10), r = fill(r, length(idx)))
                    #Push PN
                    for k in 1:length(p.l), kk in 1:length(p.q)
                        @views insertcols!(df_temp, col_str[k, kk] => round.(compute_box_iprs(ltc_p, psi[:, idx], boxidx[k], q = p.q[kk]), sigdigits = 10))
                    end
                    append!(df[x], df_temp)
                    er = false
                catch e
                    println("There was an error: ", e.msg)
                    er_num += 1
                    continue
                end
            end
        end
        CSV.write("L$(p.L)_Th$(j)_W$(jj)_E$(jjj).csv", vcat(df...))
        rm("L$(p.L)_Th$(j)_W$(jj)_E$(jjj)_temp_$(tag).csv")
    end
end

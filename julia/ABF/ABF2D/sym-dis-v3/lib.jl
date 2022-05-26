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
include("./params.jl") # read parameters from configuration file

"""
Linear map for A-IE
"""
function shift_invert_linear_map(A, E; c = 1.)
    F = factorize(c*A-E*I(size(A, 1)))
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), size(A, 1), ismutating = true, ishermitian = true)
end


function ldiv2!(y, F, x)
    y .= F\x
end

""" 
Create onsite disorder realization array in a symplectic way.
"""
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
    lists = glob("L$(p.L)/L*_Th*_W*_E*.csv")
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

function generate_fn(L, j, jj, jjj)
    return "L$(L)_Th$(j)_W$(jj)_E$(jjj)" 
end

function get_scan_idx(p::Params, scan_id::UInt64)
    # Search for unscanned index
    @label Search

    j, jj, jjj = missing_idx_finder(p)
    if (j, jj, jjj) == (0, 0, 0)
        return (j, jj, jjj)	
    end

    fn = "L$(p.L)/"*generate_fn(p.L,j, jj, jjj)
    CSV.write(fn*"_temp_$(scan_id).csv", DataFrame())
    scan_overlap = glob(fn*"_temp_*.csv")
    sleep(3rand())
    if length(scan_overlap) != 1
        println("Overlap detected")
        rm(fn*"_temp_$(scan_id).csv")
        @goto Search
    end

    return j, jj, jjj
end

struct EParams
    E_bw::Array{Float64}
    E_c::Array{Float64}
    E_del::Array{Float64}
end

function abf2d_scan(p::Params, p_E::EParams)
    mkpath("L$(p.L)")
    col_str = generate_col_names(p)
    nt = Threads.nthreads()
    rng = [MersenneTwister(p.seed + i) for i in 1:nt]
    ltc = Lattice2D(p.L, p.L, 4)
    ltc_p = Lattice2D(p.L, p.L, 2)
    boxidx = [box_inds(ltc_p, li) for li in p.l]
    scan_id = rand(UInt64)
    scan_flag = true
    while scan_flag
        j, jj, jjj = get_scan_idx(p, scan_id)
        if (j, jj, jjj) == (0, 0, 0)
            println("Nothing to scan. Terminate")
            break
        end

        #---- Energy parameters ----#
        E_c = p_E.E_c[j, jj, jjj]
        E_del = p_E.E_del[j, jj]

        #----------------------------- START DIAGONALIZATION -----------------------------#
        println("start diagonalizing at index: ", (j, jj, jjj))
        H, U = ham_fe(ltc, -2., 0., p.θ[j]) # Fully entangled hamiltonian
        

        #---- Allocate Empty dataframe ----#
        df = [DataFrame(E = Float64[], r = Int64[]) for i in 1:nt]
        for t in 1:nt, k in 1:length(p.l), kk in 1:length(p.q)
            insertcols!(df[t], col_str[k, kk] => Float64[])
        end

        @Threads.threads for r in 1:p.R÷p.nev # Realizations
            er = true
            er_num = 0
            while er
                er_num == 5 && error("5 attempts faild.")
                try
                    x = Threads.threadid()
                    
                    #---- Create scale free model ----#
                    H_dis = makesym2d(ltc, H, p.V1, p.V2, rng = rng[x])
                    D = Diagonal(dis(p.L^2, p.W[jj], rng[x]))
                    @views H_prj = project(U'*(H_dis + D)*U)
                    droptol!(H_prj, 1E-14)
                    H_prj = Hermitian(H_prj)

                    #---- Lanczos method with shift-and-invert method ----#
                    lmap = shift_invert_linear_map(H_prj, E_c) 
                    e_inv, psi, info = eigsolve(lmap, 2p.L^2, p.nev, :LM, ishermitian = true, krylovdim = max(30, 2p.nev + 1));
                    e = 1 ./ (p.L^2*real.(e_inv)) .+ E_c
                    psi = reduce(hcat, psi)

                    #---- Crop energies outside the energy bins ----#
                    idx = findall(x -> (E_c - E_del) < x < (E_c + E_del), e)
                    @views df_temp = DataFrame(E = round.(e[idx], sigdigits = 10), r = fill(r, length(idx)))

                    # Push PN
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
        CSV.write("L$(p.L)/"*generate_fn(p.L, j, jj, jjj)*".csv", vcat(df...))
        rm("L$(p.L)/"*generate_fn(p.L, j, jj, jjj)*"_temp_$(scan_id).csv")
    end
end

function abf2d_scan(p::Params)
    mkpath("L$(p.L)")
    col_str = generate_col_names(p)
    nt = Threads.nthreads()
    rng = [MersenneTwister(p.seed + i) for i in 1:nt]
    ltc = Lattice2D(p.L, p.L, 4)
    ltc_p = Lattice2D(p.L, p.L, 2)
    boxidx = [box_inds(ltc_p, li) for li in p.l]
    scan_id = rand(UInt64)
    scan_flag = true
    E_c = p.E 
    if length(E_c) > 1
        E_del = (E_c[2] - E_c[1])/p.E_bin_width
    else
        E_del = 1/p.E_bin_width 
    end
    while scan_flag
        j, jj, jjj = get_scan_idx(p, scan_id)
        if (j, jj, jjj) == (0, 0, 0)
            println("Nothing to scan. Terminate")
            break
        end

        #----------------------------- START DIAGONALIZATION -----------------------------#
        println("start diagonalizing at index: ", (j, jj, jjj))
        H, U = ham_fe(ltc, -2., 0., p.θ[j]) # Fully entangled hamiltonian
        

        # Allocate Empty dataframe
        df = [DataFrame(E = Float64[], r = Int64[]) for i in 1:nt]
        for t in 1:nt, k in 1:length(p.l), kk in 1:length(p.q)
            insertcols!(df[t], col_str[k, kk] => Float64[])
        end

        @Threads.threads for r in 1:p.R÷p.nev # Realizations
            er = true
            er_num = 0
            while er
                er_num == 5 && error("5 attempts faild.")
                try
                    x = Threads.threadid()
                    
                    #---- Create scale free model ----#
                    H_dis = makesym2d(ltc, H, p.V1, p.V2, rng = rng[x])
                    D = Diagonal(dis(p.L^2, p.W[jj], rng[x]))
                    @views H_prj = project(U'*(H_dis + D)*U)
                    droptol!(H_prj, 1E-14)
                    H_prj = Hermitian(H_prj)

                    #---- Lanczos method with shift-and-invert method ----#
                    lmap = shift_invert_linear_map(H_prj, E_c[jjj], 2p.L^2) 
                    e_inv, psi, info = eigsolve(lmap, 2p.L^2, p.nev, :LM, ishermitian = true, krylovdim = max(30, 2p.nev + 1));
                    e = 1 ./ (2p.L^2*real.(e_inv)) .+ E_c[jjj]
                    psi = reduce(hcat, psi)

                    #---- Crop energies outside the energy bins ----#
                    idx = findall(x -> (E_c[jjj] - E_del) < x < (E_c[jjj] + E_del), e)
                    @views df_temp = DataFrame(E = round.(e[idx], sigdigits = 10), r = fill(r, length(idx)))

                    # Push PN
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

        CSV.write("L$(p.L)/"*generate_fn(p.L, j, jj, jjj)*".csv", vcat(df...))
        rm("L$(p.L)/"*generate_fn(p.L, j, jj, jjj)*"_temp_$(scan_id).csv")
    end
end

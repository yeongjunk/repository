
using LinearAlgebra, SparseArrays, Arpack, KrylovKit
using LinearMaps
using Random
using DataFrames, CSV
using ArgParse, JSON
# Custom modules
using ABFSym
using Lattice
using PN

include("./ed-sf-sym-fixed-eps-params.jl") # read parameters from configuration file

function construct_linear_map(A)
    F = factorize(A)
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), size(A, 1), ismutating = true, ishermitian = true)
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


function symp_v(V1, V2, ϕ)
    V = [V1 exp(-im*ϕ)*V2; -exp(im*ϕ)*V2 V1]
    return V
end

function makesym2d(ltc, H, V1, V2)
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
            Vx = symp_v(V1, V2, 0)
            symplectic_coupling!(ltc, t0[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m, n), i, j, Id)
            symplectic_coupling!(ltc, t1[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m, n + 1), i, j, Vx)
            symplectic_coupling!(ltc, t1[2(i-1) + (j-1) + 1], I, J, val, (m, n + 1), (m, n), j, i, Vx')
            symplectic_coupling!(ltc, t2[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m + 1, n), i, j, Id)
            symplectic_coupling!(ltc, t2[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n), (m, n), j, i, Id)
            symplectic_coupling!(ltc, t3[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m + 1, n + 1), i, j, Id)
            symplectic_coupling!(ltc, t3[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n + 1), (m, n), j, i, Id)
            symplectic_coupling!(ltc, t4[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n), (m, n + 1), i, j, Id)
            symplectic_coupling!(ltc, t4[2(i-1) + (j-1) + 1], I, J, val, (m, n + 1), (m + 1, n), j, i, Id)
        end
    end
    return sparse(I, J, val, size(H, 1), size(H, 2))
end

function ham_fe_sym(ltc, Ea, Eb, θ, V1, V2)
    H_fd = ham_fd(ltc, -2., 0.)
    U1 = LUT(ltc, 0.20)
    R1 = redef1(ltc)
    R2 = redef2(ltc)
    U_first = U1*R1*U1
    U_rest = U1*R2
    H_sd = U_first*H_fd*U_first'
    H_sd_sc = makesym2d(ltc, H_sd, V1, V2)
    H_fe_sc = U_rest*H_sd_sc*U_rest'
    return H_fe_sc, U_rest*U_first
end

function estimate_bw(p, θ, W, L, rng)
    ltc = Lattice2D(L, L, 4)
    H, U = ham_fe_sym(ltc, -2, 0, θ, p.V1, p.V2) # Fully entangled hamiltonian
    D = Diagonal(dis(L*L, W, rng))
    @views H_dis = (H + D)
    vals, psi, info = eigsolve(Hermitian(H_dis), size(H_dis, 1), 1, :LM, ishermitian = true, krylovdim = 30)
    return 2*abs(maximum(vals))
end

function estimate_energy_params(p, θ, W, L, E_crop, E_bin_width, rng)
    BW = estimate_bw(p, θ, W, L, rng)
    E_c = range(0.0001, BW/2*E_crop, length = p.E_num)
    E_del = (E_c[2] - E_c[1])/E_bin_width
    return BW, E_c, E_del
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

function abf3d_scan(p::Params)
    q_str = "q".*string.(p.q)
    nt = Threads.nthreads()
    rng = [MersenneTwister(p.seed + i) for i in 1:nt]
    ltc = Lattice2D(p.L, p.L, 4)
    ltc_p = Lattice2D(p.L, p.L, 2)
    boxidx = box_inds(ltc, p.l)
    for j in 1:length(p.θ)
        fn = "L$(p.L)_Th$(j)" #File name
        H, U = ham_fe_sym(ltc, -2, 0, p.θ[j], p.V1, p.V2) # Fully entangled hamiltonian
        for jj in 1:length(p.W)
            BW, E_c, E_del = estimate_energy_params(p, p.θ[j], p.W[jj], 50, 0.9, 8, rng[1])
            BW = min(p.W[jj], 2)
            E_c = range(0.00001, BW/2*0.9, length = p.E_num)
            E_del= (E_c[2] - E_c[1])/2
            println("Disorder strength", p.W)
            println("Number of threads: $(nt)")
            println("Bandwidth estimated: $(BW)")
            println("Energy bin width: $(E_del)")
            for jjj in p.start_E_ind:p.end_E_ind
                df = [DataFrame(E = Float64[], r = Int64[]) for i in 1:nt]
                for t in 1:nt, k in 1:length(p.q)
                    insertcols!(df[t], q_str[k] => Float64[])
                end
                @Threads.threads for r in 1:p.R÷p.nev # Realizations
                    er = true
                    er_num = 0
                    while er
                        if er_num == 5; error("3 attempts faild."); end
                        try
                            x = Threads.threadid()
                            H_dis =  H + Diagonal(dis(p.L*p.L, p.W[jj], rng[x]))
                            droptol!(H_dis, 1E-12)
                            e_inv, psi, info = eigsolve(construct_linear_map(Hermitian(H_dis - E_c[jjj]*I(size(H_dis, 1)))), size(H_dis, 1),
                                p.nev, :LM, ishermitian = true, krylovdim = max(30, 2p.nev + 1));

                            e = 1 ./ real.(e_inv) .+ E_c[jjj]
                            psi = reduce(hcat, psi)
                            #Crop energies that are outside the energy bins
                            idx = findall(x -> (E_c[jjj] - E_del) <  x < (E_c[jjj] + E_del), e)
                            @views df_temp = DataFrame(E = round.(e, sigdigits = 12), r = fill(r, length(e)))
                            #Push PN
                            for k in 1:length(p.q)
                                @views insertcols!(df_temp, q_str[k] => round.(compute_box_iprs(ltc_p, psi, boxidx, q = p.q[k]), sigdigits = 12))
                            end
                            append!(df[x], df_temp)
                            er = false
                        catch e
                            println("There was an error, ", e.msg)
                            er_num += 1
                            continue
                        end
                    end
                end
                CSV.write(fn*"_W$(jj)_E$(jjj).csv", vcat(df...))
            end
        end
    end

end

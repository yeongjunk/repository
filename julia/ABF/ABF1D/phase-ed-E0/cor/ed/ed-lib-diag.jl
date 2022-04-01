using LinearAlgebra, SparseArrays
using LinearMaps
using Random
# Custom modules
using ABF
using Lattice
using PN
include("./ed-params.jl") # read parameters from configuration file
function findnearest(v, t)
    return findmin(abs.(v.-t))
end
function phase_dis(H; V::Float64 = 1., rng = nothing)
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

function eig_corr(x, i, len_r)
    N = length(x)
    y = zeros(eltype(x), len_r)
    for r in 1:len_r
        y[r] += abs.(x[i] .- x[i + r])
    end
    return y
end

function abf_scan(p::Params)
    rng = MersenneTwister(p.seed)
    ltc = Lattice1D(p.L, 2)
    ltc_p = Lattice1D(p.L,  1)
    begin_cor = p.L÷4
    end_cor = p.L÷4*3
    len_r = p.L - end_cor 
    logpsi = zeros(Float64, len_r, length(p.θ))
    logpsi_sq = zeros(Float64, len_r, length(p.θ))
    @Threads.threads for j in 1:length(p.θ) 
        H, U = ham_fe(ltc, -2, 0, p.θ[j]) # Fully entangled hamiltonian
        H = convert.(ComplexF64, H)
        logpsi_j = zeros(Float64, len_r)
        logpsi_sq_j = zeros(Float64, len_r)
        for r in 1:p.R # Realizations
            D_phase = phase_dis(H, V = 1., rng = rng)
            @views H_prj = project(U'*(p.V*D_phase)*U)
            if !p.periodic
                H_prj[1:4, end-4:end] .= 0
                H_prj[end-4:end, 1:4] .= 0
            end 
            e, vecs = eigen!(Hermitian(Matrix(H_prj)), -p.E_window, p.E_window)
            e0val, e0idx = findnearest(e, 0.)
            @views logpsi_j_r = log.(abs.(vecs[:, e0idx]))
            for i in begin_cor:end_cor
                logpsi_j .+= eig_corr(logpsi_j_r, i, len_r)
                logpsi_sq_j .+= eig_corr(logpsi_j_r.^2, i, len_r)
            end
        end
        logpsi[:, j] .= logpsi_j            
        logpsi_sq[:, j] .= logpsi_sq_j
    end
    logpsi ./= (p.R*len_r)
    logpsi_sq ./= (p.R*len_r)
    logpsi_std = sqrt.(logpsi_sq .- logpsi.^2)
    logpsi_ste = logpsi_std / sqrt(p.R*len_r)
    log10_logpsi_ste = logpsi_ste ./ (logpsi*log(10))
    fn = "output.mat"
    description = "Correlation of wavefunction at E = 0, ABF with phase disorder. First axis is the site, second axis is the theta"
    return Dict("R" => p.R, "L" => p.L, "th" => collect(p.θ),  "description" => description, "logpsi_mean" => logpsi, "logpsi_std" => logpsi_std, "logpsi_ste" => logpsi_ste, "log10_logpsi_ste" => log10_logpsi_ste)
end

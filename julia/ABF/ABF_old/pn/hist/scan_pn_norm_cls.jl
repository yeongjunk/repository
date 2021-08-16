using MAT
using Dates
using JSON
using ArgParse
using StaticArrays

include("library/abf2_pnscan.jl")
include("library/abf2_cob.jl")

@doc """
Calculate pn, norm of psi in cls.
"""
function cls_scan(p::Parameters)
    W, θ, _ = expand_params(p)
    pn = zeros(Float64, p.θI_num, p.W_num)
    norm_a = Array{Float64}(undef,p.θI_num, p.W_num)
    norm_b = Array{Float64}(undef,p.θI_num, p.W_num)

    counters = zeros(Int64, p.θI_num, p.W_num) # of samples for each W, θ
    nt = Threads.nthreads()
    H = zeros(Float64, 2p.N, 2p.N, nt) # individual hamiltonian matrix storage for each threads.

    for n in 1:p.R #Iteration for each disorder realizations
        rng = MersenneTwister(p.seed+n)
        r_arr = rand(rng, 2p.N) .- 0.5
        for i in 1:length(θ)
            for t = 1:nt
                overwrite_ham_abf2!(view(H,:, :, t), θ[i], θ[i], p.N)
            end
                Threads.@threads for j in 1:length(W) # Iteration for W
                    H_dis = view(H, :,:,Threads.threadid())
                    r_onsite = W[j]*r_arr
                    add_diag!(H_dis, r_onsite)
                    eig = eigen(Hermitian(H_dis));
                    add_diag!(H_dis, -r_onsite) #remove disorder for recycling


                    indices = findinrange(eig.values, p.E, 2p.E_range)
                    n_indices = length(indices)

                    if n_indices == 0 # there is no sample within energy window
                        println("no state, i=$(i), j=$(j). Skips iteration")
                        continue
                    end

                    counters[i,j] += n_indices # number of samples collected within energy window.

                    t_pn = 0
                    t_norm_a = 0
                    t_norm_b = 0
                    U = unitary(-θ[i])
                    for idx in indices
                         a, b = splitpsi(eig.vectors[:,idx])
                         a, b = unitary(a, b, U)

                         a, b = uc_redef_n(a, b)

                         a, b = unitary(a, b, U)
                         t_norm_a += norm(a)
                         t_norm_b += norm(b)
                         t_pn += compute_pn(abs.(vcat(a,b)))
                    end
                    pn[i,j] += t_pn
                    norm_a[i, j] +=t_norm_a
                    norm_b[i, j] +=t_norm_b
            end

        end
    end
    pn = pn./counters
    norm_a = norm_a./counters
    norm_b = norm_b./counters

    println("mincount = $(findmin(counters)), maxcount = $(findmax(counters))")
    return pn, norm_a, norm_b
end



function main(args)
    opts = ArgParseSettings(description="Scan and compute pn for all parameters of nu=2 ABF")
    @add_arg_table! opts begin
    "c"
        help = "configuration"
        #nargs = 1
        arg_type = AbstractString
    end
    # Parse the arguments
    popts   = parse_args(opts)     # parse the arguments
    config  = JSON.parsefile(popts["c"])
    params = readconfig(config)   # create the Parameters struct from Dict config
    BLAS.set_num_threads(1) #Turn off BLAS multi-threading
    @time pn, norm_a, norm_b = cls_scan(params)
    config["pn"]  = pn
    config["norm_a"] = norm_a
    config["norm_b"] = norm_b

    datestring =Dates.format(now(), "ddmmyy-HHMMSS") # Datetime as output filename
    fn_out = datestring*".mat"
    matwrite(fn_out, config) # write compressed MATLAB file
end

main(ARGS)

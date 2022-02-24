using MAT
using Dates
@everywhere include("/home/yjkim/codes/repository/julia/module/NNQuasi1D/src/tmm.jl")
@everywhere using Random
@everywhere using LinearAlgebra
@everywhere using SparseArrays
@everywhere using Distributed
@everywhere LinearAlgebra.BLAS.set_num_threads(1)
@everywhere W = [10., 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0,
    14.3, 14.5, 15.0, 15.3, 15.5, 15.8, 16.0, 16.1, 16.2, 16.3, 16.4,
    16.5, 16.6, 16.7, 16.8, 17.0, 17.2, 17.5, 18.0, 18.5, 19.0,
    20.0, 21.0, 22.0, 24.0, 26.0, 28.0, 30.0]
@everywhere M = collect(3:1:10)
@everywhere N = Int64(10^7)
@everywhere seed = [9951 + i for i in 1:100]
for i in 1:nprocs()
	@spawnat procs()[i] rng = MersenneTwister(seed[i])
end

println("N = ", N)
println("Number of Processors: ",nprocs())

@time xi = pmap(t -> ((w, m)=t; tmm_bar(M = m, E = 0., W = w, N = N, N_qr = 10)), Iterators.product(W, M))

#for i in 1:length(M)
#    xi[:, i] /= M[i]
#end
display(xi)
MAT.matwrite(Dates.format(now(), "HHMMSS")*"-xi-cubic.mat", Dict("xi" => xi, "W" => W, "M" => M, "N"=>N, "seed"=>seed, "N_qr" => 10))

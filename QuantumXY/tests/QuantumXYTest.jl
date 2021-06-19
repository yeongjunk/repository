module QuantumXYTest

push!(LOAD_PATH, homedir()*"/codes/project/ising/QuantumXY/src/")
using QuantumXY
using SparseArrays
using LinearAlgebra
using Test
using Random

Random.seed!(1234321)
up = [1.; 0]
down = [0.; 1]

@testset "Quantum XY model" verbose = true begin
    @testset "Kronecker Product" begin
        N = 5; J = rand(); h = rand()
        H_kron = ham_xy_kron(N = N, J = J, h = h)
        H = ham_xy(N = N, J = J, h = h)
        @test isequal(H, H_kron)
    end

    @testset "GS for J = 0, h > 0 " begin
        H = ham_xy(N = 4, J = 0, h = 1.)
        values, vectors = eigen(Hermitian(Matrix(H)))
        @test vectors[:, 1] == kron(up, up, up, up)
    end
    @testset "GS for J = 0, h < 0 " begin
        H = ham_xy(N = 4, J = 0, h = -1.)
        values, vectors = eigen(Hermitian(Matrix(H)))
        @test vectors[:, 1] == kron(down, down, down, down)
    end
end

end #module

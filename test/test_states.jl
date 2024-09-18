@testitem "States" begin
    using Gabs
    using StaticArrays

    @testset "vacuum states" begin
        @test vacuumstate() isa GaussianState
        @test vacuumstate(SVector{2}, SMatrix{2,2}) isa GaussianState
    end

    @testset "thermal states" begin
        n = rand(Int64)
        @test thermalstate(n) isa GaussianState
        @test thermalstate(SVector{2}, SMatrix{2,2}, n) isa GaussianState
    end

    @testset "coherent states" begin
        alpha = rand(ComplexF64)
        @test coherentstate(alpha) isa GaussianState
        @test coherentstate(SVector{2}, SMatrix{2,2}, alpha) isa GaussianState
    end

    @testset "squeezed states" begin
        r, theta = rand(Float64), rand(Float64)
        @test squeezedstate(r, theta) isa GaussianState
        @test squeezedstate(SVector{2}, SMatrix{2,2}, r, theta) isa GaussianState
    end

    @testset "epr states" begin
        r, theta = rand(Float64), rand(Float64)
        @test eprstate(r, theta) isa GaussianState
        @test eprstate(SVector{4}, SMatrix{4,4}, r, theta) isa GaussianState
    end

    @testset "tensor products" begin
        v = vacuumstate()
        vs = tensor(v, v)
        @test vs isa GaussianState
        @test tensor(SVector{4}, SMatrix{4,4}, v, v) isa GaussianState
        @test vs == v ⊗ v

        alpha = rand(ComplexF64)
        c = coherentstate(alpha)
        @test tensor(c, tensor(v, v)) == c ⊗ v ⊗ v


        vstatic = vacuumstate(SVector{2}, SMatrix{2,2})
        tpstatic = vstatic ⊗ vstatic ⊗ vstatic
        @test tpstatic.mean isa SVector{6}
        @test tpstatic.covar isa SMatrix{6,6}
        tp = vstatic ⊗ v ⊗ vstatic
        @test tp.mean isa Vector
        @test tp.covar isa Matrix
    end

    @testset "partial trace" begin
        alpha = rand(Float64)
        r, theta = rand(Float64), rand(Float64)
        n = rand(Int)
        s1, s2, s3 = coherentstate(alpha), squeezedstate(r, theta), thermalstate(n)
        state = s1 ⊗ s2 ⊗ s3
        @test ptrace(state, 1) == s1
        @test ptrace(state, 2) == s2
        @test ptrace(state, 3) == s3
        @test ptrace(state, [1, 2]) == s1 ⊗ s2
        @test ptrace(state, [1, 3]) == s1 ⊗ s3
        @test ptrace(state, [2, 3]) == s2 ⊗ s3

        @test ptrace(SVector{2}, SMatrix{2,2}, state, 1) isa GaussianState
        @test ptrace(SVector{4}, SMatrix{4,4}, state, [1, 3]) isa GaussianState
end
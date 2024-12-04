@testitem "Quadrature pair basis" begin
    using Gabs
    using StaticArrays

    basis1 = QuadPairBasis(1)

    @testset "vacuum states" begin
        @test vacuumstate(basis1) isa GaussianState
        @test vacuumstate(SVector{2}, SMatrix{2,2}, basis1) isa GaussianState
    end

    @testset "thermal states" begin
        n = rand(Int64)
        @test thermalstate(basis1, n) isa GaussianState
        @test thermalstate(SVector{2}, SMatrix{2,2}, basis1, n) isa GaussianState
    end

    @testset "coherent states" begin
        alpha = rand(ComplexF64)
        @test coherentstate(basis1, alpha) isa GaussianState
        @test coherentstate(SVector{2}, SMatrix{2,2}, basis1, alpha) isa GaussianState
    end

    @testset "squeezed states" begin
        r, theta = rand(Float64), rand(Float64)
        @test squeezedstate(basis1, r, theta) isa GaussianState
        @test squeezedstate(SVector{2}, SMatrix{2,2}, basis1, r, theta) isa GaussianState
    end

    @testset "epr states" begin
        r, theta = rand(Float64), rand(Float64)
        @test eprstate(2*basis1, r, theta) isa GaussianState
        @test eprstate(SVector{4}, SMatrix{4,4}, 2*basis1, r, theta) isa GaussianState
    end

    @testset "tensor products" begin
        v = vacuumstate(basis1)
        vs = tensor(v, v)
        @test vs isa GaussianState
        @test tensor(SVector{4}, SMatrix{4,4}, v, v) isa GaussianState
        @test vs == v ⊗ v
        @test isapprox(vs, v ⊗ v)

        alpha = rand(ComplexF64)
        c = coherentstate(basis1, alpha)
        @test tensor(c, tensor(v, v)) == c ⊗ v ⊗ v

        vstatic = vacuumstate(SVector{2}, SMatrix{2,2}, basis1)
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
        s1, s2, s3 = coherentstate(basis1, alpha), squeezedstate(basis1, r, theta), thermalstate(basis1, n)
        state = s1 ⊗ s2 ⊗ s3
        @test ptrace(state, 1) == s1
        @test ptrace(state, 2) == s2
        @test ptrace(state, 3) == s3
        @test ptrace(state, [1, 2]) == s1 ⊗ s2
        @test ptrace(state, [1, 3]) == s1 ⊗ s3
        @test ptrace(state, [2, 3]) == s2 ⊗ s3

        sstatic = coherentstate(SVector{2}, SMatrix{2,2}, basis1, alpha)
        tpstatic = sstatic ⊗ sstatic ⊗ sstatic
        @test ptrace(tpstatic, 1) == sstatic
        @test ptrace(tpstatic, [1,3]) == sstatic ⊗ sstatic

        @test ptrace(SVector{2}, SMatrix{2,2}, state, 1) isa GaussianState
        @test ptrace(SVector{4}, SMatrix{4,4}, state, [1, 3]) isa GaussianState
    end
end
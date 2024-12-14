@testitem "States" begin
    import Gabs: _changebasis
    using Gabs
    using StaticArrays

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)

    @testset "vacuum states" begin
        state = vacuumstate(qpairbasis)
        @test state isa GaussianState
        @test vacuumstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis) isa GaussianState
        @test vacuumstate(qblockbasis) == _changebasis(state, QuadBlockBasis)
    end

    @testset "thermal states" begin
        n = rand(Int64)
        state = thermalstate(qpairbasis, n)
        @test state isa GaussianState
        @test thermalstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, n) isa GaussianState
        @test thermalstate(qblockbasis, n) == _changebasis(state, QuadBlockBasis)
    end

    @testset "coherent states" begin
        alpha = rand(ComplexF64)
        state = coherentstate(qpairbasis, alpha)
        @test state isa GaussianState
        @test coherentstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, alpha) isa GaussianState
        @test coherentstate(qblockbasis, alpha) == _changebasis(state, QuadBlockBasis)
    end

    @testset "squeezed states" begin
        r, theta = rand(Float64), rand(Float64)
        state = squeezedstate(qpairbasis, r, theta)
        @test state isa GaussianState
        @test squeezedstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, r, theta) isa GaussianState
        @test squeezedstate(qblockbasis, r, theta) == _changebasis(state, QuadBlockBasis)
    end

    @testset "epr states" begin
        r, theta = rand(Float64), rand(Float64)
        state = eprstate(2*qpairbasis, r, theta)
        @test state isa GaussianState
        @test eprstate(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, r, theta) isa GaussianState
        @test eprstate(2*qblockbasis, r, theta) == _changebasis(state, QuadBlockBasis)
    end

    @testset "tensor products" begin
        v = vacuumstate(qpairbasis)
        vs = tensor(v, v)
        @test vs isa GaussianState
        @test tensor(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, v, v) isa GaussianState
        @test vs == v ⊗ v
        @test isapprox(vs, v ⊗ v)

        alpha = rand(ComplexF64)
        c = coherentstate(qpairbasis, alpha)
        @test tensor(c, tensor(v, v)) == c ⊗ v ⊗ v

        vstatic = vacuumstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis)
        tpstatic = vstatic ⊗ vstatic ⊗ vstatic
        @test tpstatic.mean isa SVector{6*nmodes}
        @test tpstatic.covar isa SMatrix{6*nmodes,6*nmodes}
        tp = vstatic ⊗ v ⊗ vstatic
        @test tp.mean isa Vector
        @test tp.covar isa Matrix
    end

    @testset "partial trace" begin
        qpairbasis1 = QuadPairBasis(1)
        alpha = rand(Float64)
        r, theta = rand(Float64), rand(Float64)
        n = rand(Int)
        s1, s2, s3 = coherentstate(qpairbasis1, alpha), squeezedstate(qpairbasis1, r, theta), thermalstate(qpairbasis1, n)
        state = s1 ⊗ s2 ⊗ s3
        @test ptrace(state, 1) == s1
        @test ptrace(state, 2) == s2
        @test ptrace(state, 3) == s3
        @test ptrace(state, [1, 2]) == s1 ⊗ s2
        @test ptrace(state, [1, 3]) == s1 ⊗ s3
        @test ptrace(state, [2, 3]) == s2 ⊗ s3

        sstatic = coherentstate(SVector{2}, SMatrix{2,2}, qpairbasis1, alpha)
        tpstatic = sstatic ⊗ sstatic ⊗ sstatic
        @test ptrace(tpstatic, 1) == sstatic
        @test ptrace(tpstatic, [1,3]) == sstatic ⊗ sstatic

        @test ptrace(SVector{2}, SMatrix{2,2}, state, 1) isa GaussianState
        @test ptrace(SVector{4}, SMatrix{4,4}, state, [1, 3]) isa GaussianState
    end
end
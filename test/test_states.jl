@testitem "States" begin
    using Gabs
    using StaticArrays
    using LinearAlgebra: det

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)

    @testset "vacuum states" begin
        state, array_state, static_array, static_state = vacuumstate(qpairbasis), vacuumstate(Array, qpairbasis),vacuumstate(SArray, qpairbasis), vacuumstate(SVector, SMatrix, qpairbasis)
        @test state isa GaussianState && array_state isa GaussianState && static_state isa GaussianState
        @test vacuumstate(qblockbasis) == changebasis(QuadBlockBasis, state)
        @test state.ħ == 2 && array_state.ħ == 2 && static_state.ħ == 2
    end

    @testset "thermal states" begin
        n = rand(1:5)
        ns = rand(1:5, nmodes)
        state_pair = thermalstate(qpairbasis, n)
        state_block = thermalstate(qblockbasis, n)
        @test state_pair isa GaussianState && state_block isa GaussianState
        @test thermalstate(Array, qpairbasis, n) isa GaussianState
        @test thermalstate(SArray, qpairbasis, n) isa GaussianState
        @test thermalstate(SVector, SMatrix, qpairbasis, n) isa GaussianState
        @test thermalstate(qblockbasis, n) == changebasis(QuadBlockBasis, state_pair)
        @test state_pair == changebasis(QuadPairBasis, state_block) && state_block == changebasis(QuadBlockBasis, state_pair)
        @test state_pair == changebasis(QuadPairBasis, state_pair) && state_block == changebasis(QuadBlockBasis, state_block)
        @test thermalstate(qblockbasis, ns) == changebasis(QuadBlockBasis, thermalstate(qpairbasis, ns))
        @test isgaussian(state_pair, atol = 1e-4)
        @test state_pair.ħ == 2 && state_block.ħ == 2
    end

    @testset "coherent states" begin
        alpha = rand(ComplexF64)
        alphas = rand(ComplexF64, nmodes)
        state_pair = coherentstate(qpairbasis, alpha)
        state_block = coherentstate(qblockbasis, alpha)
        @test state_pair isa GaussianState && state_block isa GaussianState
        @test coherentstate(Array, qpairbasis, alpha) isa GaussianState
        @test coherentstate(SArray, qpairbasis, alpha) isa GaussianState
        @test coherentstate(SVector, SMatrix, qpairbasis, alpha) isa GaussianState
        @test coherentstate(qblockbasis, alpha) == changebasis(QuadBlockBasis, state_pair)
        @test state_pair == changebasis(QuadPairBasis, state_block) && state_block == changebasis(QuadBlockBasis, state_pair)
        @test state_pair == changebasis(QuadPairBasis, state_pair) && state_block == changebasis(QuadBlockBasis, state_block)
        @test coherentstate(qblockbasis, alphas) == changebasis(QuadBlockBasis, coherentstate(qpairbasis, alphas))
        @test isgaussian(state_pair, atol = 1e-4)
        @test state_pair.ħ == 2 && state_block.ħ == 2
    end

    @testset "squeezed states" begin
        r, theta = rand(Float64), rand(Float64)
        rs, thetas = rand(Float64, nmodes), rand(Float64, nmodes)
        state, array_state, static_array, static_state = squeezedstate(qpairbasis, r, theta), squeezedstate(Array, qpairbasis, r, theta),squeezedstate(SArray, qpairbasis, r, theta),  squeezedstate(SVector, SMatrix, qpairbasis, r, theta)
        @test state isa GaussianState && array_state isa GaussianState && static_state isa GaussianState
        @test squeezedstate(qblockbasis, r, theta) == changebasis(QuadBlockBasis, state)
        @test squeezedstate(qblockbasis, rs, thetas) == changebasis(QuadBlockBasis, squeezedstate(qpairbasis, rs, thetas))
        @test state.ħ == 2 && array_state.ħ == 2 && static_state.ħ == 2
    end

    @testset "epr states" begin
        r, theta = rand(Float64), rand(Float64)
        rs, thetas = rand(Float64, nmodes), rand(Float64, nmodes)
        state, array_state, static_array, static_state = eprstate(2*qpairbasis, r, theta), eprstate(Array, 2*qpairbasis, r, theta), eprstate(SArray, 2*qpairbasis, r, theta),  eprstate(SVector, SMatrix, 2*qpairbasis, r, theta)
        @test state isa GaussianState && array_state isa GaussianState && static_state isa GaussianState
        @test eprstate(SVector, SMatrix, 2*qpairbasis, r, theta) isa GaussianState
        @test eprstate(2*qblockbasis, r, theta) == changebasis(QuadBlockBasis, state)
        @test eprstate(2*qblockbasis, rs, thetas) == changebasis(QuadBlockBasis, eprstate(2*qpairbasis, rs, thetas))
        @test state.ħ == 2 && array_state.ħ == 2 && static_state.ħ == 2
    end

    @testset "tensor products" begin
        v = vacuumstate(qpairbasis)
        vs = tensor(v, v)
        @test vs isa GaussianState
        @test tensor(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, v, v) isa GaussianState
        @test vs == v ⊗ v
        @test isapprox(vs, v ⊗ v, atol = 1e-10)

        alpha = rand(ComplexF64)
        c = coherentstate(qpairbasis, alpha)
        @test tensor(c, tensor(v, v)) == c ⊗ v ⊗ v

        r, theta = rand(Float64), rand(Float64)
        sq = squeezedstate(qblockbasis, r, theta)
        sqs = squeezedstate(2*qblockbasis, repeat([r], 2*nmodes), repeat([theta], 2*nmodes))
        @test sq ⊗ sq == sqs

        vstatic = vacuumstate(SVector, SMatrix, qpairbasis)
        tpstatic = vstatic ⊗ vstatic ⊗ vstatic
        @test tpstatic.mean isa SVector{6*nmodes}
        @test tpstatic.covar isa SMatrix{6*nmodes,6*nmodes}
        tp = vstatic ⊗ v ⊗ vstatic
        @test tp.mean isa Vector
        @test tp.covar isa SMatrix
    end

    @testset "partial trace" begin
        qpairbasis1 = QuadPairBasis(1)
        qblockbasis1 = QuadBlockBasis(1)
        alpha = rand(Float64)
        r, theta = rand(Float64), rand(Float64)
        n = rand(Int)
        s1_qpair, s2_qpair, s3_qpair = coherentstate(qpairbasis1, alpha), squeezedstate(qpairbasis1, r, theta), thermalstate(qpairbasis1, n)
        state_qpair = s1_qpair ⊗ s2_qpair ⊗ s3_qpair
        @test ptrace(state_qpair, 1) == s1_qpair
        @test ptrace(state_qpair, 2) == s2_qpair
        @test ptrace(state_qpair, 3) == s3_qpair
        @test ptrace(state_qpair, [1, 2]) == s1_qpair ⊗ s2_qpair
        @test ptrace(state_qpair, [1, 3]) == s1_qpair ⊗ s3_qpair
        @test ptrace(state_qpair, [2, 3]) == s2_qpair ⊗ s3_qpair

        s1_qblock, s2_qblock, s3_qblock = coherentstate(qblockbasis1, alpha), squeezedstate(qblockbasis1, r, theta), thermalstate(qblockbasis1, n)
        state_qblock = s1_qblock ⊗ s2_qblock ⊗ s3_qblock
        @test ptrace(state_qblock, 1) == s1_qblock
        @test ptrace(state_qblock, 2) == s2_qblock
        @test ptrace(state_qblock, 3) == s3_qblock
        @test ptrace(state_qblock, [1, 2]) == s1_qblock ⊗ s2_qblock
        @test ptrace(state_qblock, [1, 3]) == s1_qblock ⊗ s3_qblock
        @test ptrace(state_qblock, [2, 3]) == s2_qblock ⊗ s3_qblock

        sstatic = coherentstate(SVector, SMatrix, qpairbasis1, alpha)
        tpstatic = sstatic ⊗ sstatic ⊗ sstatic
        @test ptrace(tpstatic, 1) == sstatic
        @test ptrace(tpstatic, [1,3]) == sstatic ⊗ sstatic

        @test ptrace(SVector{2}, SMatrix{2,2}, state_qpair, 1) isa GaussianState
        @test ptrace(SVector{4}, SMatrix{4,4}, state_qpair, [1, 3]) isa GaussianState

        qpairbasis4 = QuadPairBasis(4)
        qblockbasis4 = QuadBlockBasis(4)

        eprstates_qpair = eprstate(qpairbasis4, r, theta)
        eprstates_qblock = eprstate(qblockbasis4, r, theta)

        @test ptrace(eprstates_qpair, [1, 2]) == eprstate(QuadPairBasis(2), r, theta)
        @test ptrace(eprstates_qblock, [1, 2]) == eprstate(QuadBlockBasis(2), r, theta)
    end

    @testset "symplectic spectrum" begin
        nmodes = rand(1:5)
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)

        s_qpair = randstate(qpairbasis)
        s_qblock = randstate(qblockbasis)

        spec_qpair = sympspectrum(s_qpair)
        spec_qblock = sympspectrum(s_qblock)

        @test all(i > 1 || isapprox(i, 1, atol=1e-5) for i in spec_qpair)
        @test all(i > 1 || isapprox(i, 1, atol=1e-5) for i in spec_qblock)

        @test isapprox(det(s_qpair.covar), prod(abs2, spec_qpair), atol=1e-3)
        @test isapprox(det(s_qblock.covar), prod(abs2, spec_qblock), atol=1e-3)
    end
end
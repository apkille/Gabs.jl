@testitem "Symbolic States" begin
    using Gabs
    using Symbolics
    using StaticArrays
    using LinearAlgebra: det

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)

    @testset "Symbolic Partial Trace of Beamsplitter-EPR State" begin
        @variables r θ τ
        b = QuadBlockBasis(2)
        st = eprstate(b, r, θ)
        op = beamsplitter(b, τ)
        newst = ptrace(op * st, 1)
        expected_mean = zeros(2)
        τ_sqrt = sqrt(τ)
        one_minus_τ_sqrt = sqrt(1 - τ)
        expected_covariance = [
            (0.5cosh(2r) * one_minus_τ_sqrt - 0.5cos(θ) * sinh(2r) * τ_sqrt) * one_minus_τ_sqrt +
            (0.5τ_sqrt * cosh(2r) - 0.5cos(θ) * sinh(2r) * one_minus_τ_sqrt) * τ_sqrt  -sinh(2r) * sin(θ) * τ_sqrt * one_minus_τ_sqrt;
            -sinh(2r) * sin(θ) * τ_sqrt * one_minus_τ_sqrt  (0.5cosh(2r) * one_minus_τ_sqrt + 0.5cos(θ) * sinh(2r) * τ_sqrt) * one_minus_τ_sqrt +
            (0.5τ_sqrt * cosh(2r) + 0.5cos(θ) * sinh(2r) * one_minus_τ_sqrt) * τ_sqrt
        ]
        @test newst isa GaussianState
        @test newst.mean == expected_mean
        @test all(isequal.(newst.covar, expected_covariance))
    end

    @testset "Symbolic EPR states" begin
        @variables r θ
        nmodes = rand(1:5)
        @variables rs[1:nmodes] thetas[1:nmodes]
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)
        state = eprstate(2 * qpairbasis, r, θ)
        @test state isa GaussianState
        @test eprstate(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, r, θ) isa GaussianState
        @test iszero(simplify(eprstate(2*qblockbasis, r, θ).covar - changebasis(QuadBlockBasis, state).covar))
        @test iszero(simplify(eprstate(2*qblockbasis, r, θ).mean - changebasis(QuadBlockBasis, state).mean))
        state_pair = eprstate(2*qpairbasis, collect(rs), collect(thetas))
        state_block = eprstate(2*qblockbasis, collect(rs), collect(thetas))
        @test iszero(simplify(state_block.covar - changebasis(QuadBlockBasis, state_pair).covar))
        @test iszero(simplify(state_block.mean - changebasis(QuadBlockBasis, state_pair).mean))
    end

    @testset "Symbolic squeezed states" begin
        @variables r theta
        @variables rs[1:nmodes] thetas[1:nmodes]
        state = squeezedstate(qpairbasis, r, theta)
        @test state isa GaussianState
        @test squeezedstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, r, theta) isa GaussianState
        @test all(isequal.(squeezedstate(qblockbasis, r, theta).covar, changebasis(QuadBlockBasis, state).covar))
        @test all(isequal.(squeezedstate(qblockbasis, r, theta).mean, changebasis(QuadBlockBasis, state).mean))
        rs_vec = collect(rs)
        thetas_vec = collect(thetas)
        state = squeezedstate(qpairbasis, rs_vec, thetas_vec)
        @test state isa GaussianState
        @test all(isequal.(squeezedstate(qblockbasis, rs_vec, thetas_vec).covar, changebasis(QuadBlockBasis, state).covar))
        @test all(isequal.(squeezedstate(qblockbasis, rs_vec, thetas_vec).mean, changebasis(QuadBlockBasis, state).mean))
    end

    @testset "Symbolic coherent states" begin
        @variables α
        @variables alphas[1:nmodes]
        state_pair = coherentstate(qpairbasis, α)
        state_block = coherentstate(qblockbasis, α)
        @test state_pair isa GaussianState && state_block isa GaussianState
        @test coherentstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, α) isa GaussianState
        @test coherentstate(qblockbasis, α).covar == changebasis(QuadBlockBasis, state_pair).covar
        @test isequal(coherentstate(qblockbasis, α).mean, changebasis(QuadBlockBasis, state_pair).mean)
    end
end

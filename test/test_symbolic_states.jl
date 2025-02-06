@testitem "Symbolic States" begin
    using Gabs
    using Symbolics
    using StaticArrays
    using LinearAlgebra: det

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)

    @testset "Symbolic Partial Trace of Beamsplitter-EPR State" begin
        @variables r θ
        @variables rs[1:nmodes] thetas[1:nmodes]
        state = eprstate(2 * qpairbasis, r, θ)
        @test state isa GaussianState
        state_array = eprstate(2 * qpairbasis, collect(rs), collect(thetas))
        @test state_array isa GaussianState
        @test iszero(simplify(eprstate(2 * qblockbasis, r, θ).covar - changebasis(QuadBlockBasis, state).covar))
        @test iszero(simplify(eprstate(2 * qblockbasis, r, θ).mean - changebasis(QuadBlockBasis, state).mean))
        state_pair = eprstate(2 * qpairbasis, collect(rs), collect(thetas))
        state_block = eprstate(2 * qblockbasis, collect(rs), collect(thetas))
        @test iszero(simplify(state_block.covar - changebasis(QuadBlockBasis, state_pair).covar))
        @test iszero(simplify(state_block.mean - changebasis(QuadBlockBasis, state_pair).mean))
        @test eprstate(2 * qpairbasis, collect(rs), collect(thetas)) isa GaussianState
        @test all(isequal.(eprstate(2 * qblockbasis, collect(rs), collect(thetas)).covar, changebasis(QuadBlockBasis, state_pair).covar))
        @test all(isequal.(eprstate(2 * qblockbasis, collect(rs), collect(thetas)).mean, changebasis(QuadBlockBasis, state_pair).mean))
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
        alphas_vec = vcat([real(alphas[i]) for i in 1:nmodes], [imag(alphas[i]) for i in 1:nmodes])
        @test all.(isequal(coherentstate(qblockbasis, alphas_vec).mean, changebasis(QuadBlockBasis, coherentstate(qpairbasis, alphas_vec)).mean))
        @test all.(isequal(coherentstate(qblockbasis, alphas_vec).covar, changebasis(QuadBlockBasis, coherentstate(qpairbasis, alphas_vec)).covar))
    end
end

@testitem "Factorizations" begin
    using Gabs
    using LinearAlgebra

    @testset "williamson" begin
        nmodes = rand(1:5)
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)
        
        state_pair = randstate(qpairbasis)
        state_block = randstate(qblockbasis)

        F_pair = williamson(state_pair)
        S_pair, spectrum_pair = F_pair
        F_block = williamson(state_block)
        S_block, spectrum_block = F_block

        @test S_pair == F_pair.S && spectrum_pair == F_pair.spectrum
        @test S_block == F_block.S && spectrum_block == F_block.spectrum

        @test isapprox(S_pair * state_pair.covar * S_pair', Diagonal(repeat(spectrum_pair, inner = 2)), atol = 1e-5)
        @test isapprox(S_block * state_block.covar * S_block', Diagonal(repeat(spectrum_block, 2)), atol = 1e-5)

        @test issymplectic(qpairbasis, S_pair, atol = 1e-5)
        @test issymplectic(qblockbasis, S_block, atol = 1e-5)
    end
end
@testitem "Symbolic Unitaries" begin
    using Gabs
    using Symbolics
    using StaticArrays
    using LinearAlgebra: det

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)

    @testset "displacement operator" begin
        @variables α alphas[1:5]
        op_pair = displace(qpairbasis, α)
        op_block = displace(qblockbasis, α)
        @test op_pair isa GaussianUnitary && op_block isa GaussianUnitary
        @test displace(Array, qpairbasis, α) isa GaussianUnitary
        @test displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, α) isa GaussianUnitary
        @test all.(isequal(op_pair.symplectic, changebasis(QuadBlockBasis, op_pair).symplectic))
        @test all.(isequal(op_pair.symplectic, changebasis(QuadBlockBasis, op_block).symplectic))
        @test all.(isequal(displace(qblockbasis, α).disp, changebasis(QuadBlockBasis, op_pair).disp))
        @test all.(isequal(displace(qblockbasis, α).symplectic, changebasis(QuadBlockBasis, op_pair).symplectic))
        @test issymplectic(qpairbasis, op_pair.symplectic, atol = 1e-4)
        @test isgaussian(op_pair, atol = 1e-4)
        alphas_vec = vcat([real(alphas[i]) for i in 1:nmodes], [imag(alphas[i]) for i in 1:nmodes])
        op_pair = displace(qpairbasis, alphas_vec)
        op_block = displace(qblockbasis, alphas_vec)
        @test op_pair isa GaussianUnitary && op_block isa GaussianUnitary
        @test displace(Array, qpairbasis, alphas_vec) isa GaussianUnitary
        @test displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, alphas_vec) isa GaussianUnitary
        @test all.(isequal(op_pair.symplectic, changebasis(QuadBlockBasis, op_pair).symplectic))
        @test all.(isequal(op_pair.symplectic, changebasis(QuadBlockBasis, op_block).symplectic))
        @test all.(isequal(displace(qblockbasis, alphas_vec).disp, changebasis(QuadBlockBasis, op_pair).disp))
        @test all.(isequal(displace(qblockbasis, alphas_vec).symplectic, changebasis(QuadBlockBasis, op_pair).symplectic))
        @test issymplectic(qpairbasis, op_pair.symplectic, atol = 1e-4)
        @test isgaussian(op_pair, atol = 1e-4)
    end

    @testset "Symbolic squeeze operator" begin
        @variables r theta
        op_pair = squeeze(qpairbasis, r, theta)
        op_block = squeeze(qblockbasis, r, theta)
        @test op_pair isa GaussianUnitary && op_block isa GaussianUnitary
        @test squeeze(Array, qpairbasis, r, theta) isa GaussianUnitary
        @test squeeze(SVector{2*nmodes}, SMatrix{2*nmodes, 2*nmodes}, qpairbasis, r, theta) isa GaussianUnitary
        @test all.(isequal(squeeze(qblockbasis, r, theta).disp, changebasis(QuadBlockBasis, op_pair).disp))
        @test all.(isequal(squeeze(qblockbasis, r, theta).symplectic, changebasis(QuadBlockBasis, op_pair).symplectic))
        @variables rs[1:nmodes] thetas[1:nmodes]
        rs_vec = collect(rs)
        thetas_vec = collect(thetas)
        op_pair_arr = squeeze(qpairbasis, rs_vec, thetas_vec)
        op_block_arr = squeeze(qblockbasis, rs_vec, thetas_vec)
        @test op_pair_arr isa GaussianUnitary && op_block_arr isa GaussianUnitary
        @test squeeze(Array, qpairbasis, rs_vec, thetas_vec) isa GaussianUnitary
        @test squeeze(SVector{2*nmodes}, SMatrix{2*nmodes, 2*nmodes}, qpairbasis, rs_vec, thetas_vec) isa GaussianUnitary
        @test all.(isequal(squeeze(qblockbasis, rs_vec, thetas_vec).disp, changebasis(QuadBlockBasis, op_pair_arr).disp))
        @test all.(isequal(squeeze(qblockbasis, rs_vec, thetas_vec).symplectic, changebasis(QuadBlockBasis, op_pair_arr).symplectic))
    end
end

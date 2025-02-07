@testitem "Symbolic Unitaries" begin
    using Gabs
    using Symbolics
    using StaticArrays
    using LinearAlgebra: det

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)

    @testset "Symbolic displacement operator" begin
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


    for (name, op_func, factor) in [("squeeze", squeeze, 1), ("two-mode squeeze", twosqueeze, 2)]
        @testset "Symbolic $name operator" begin
            @variables r theta
            op_pair = op_func(factor * qpairbasis, r, theta)
            op_block = op_func(factor * qblockbasis, r, theta)
            @test op_pair isa GaussianUnitary && op_block isa GaussianUnitary
            @test op_func(Array, factor * qpairbasis, r, theta) isa GaussianUnitary
            @test op_func(SVector{factor * 2 * nmodes}, SMatrix{factor * 2 * nmodes, factor * 2 * nmodes}, factor * qpairbasis, r, theta) isa GaussianUnitary
            @test all.(isequal(op_func(factor * qblockbasis, r, theta).disp, changebasis(QuadBlockBasis, op_pair).disp))
            @test all.(isequal(op_func(factor * qblockbasis, r, theta).symplectic, changebasis(QuadBlockBasis, op_pair).symplectic))
            @variables rs[1:nmodes] thetas[1:nmodes]
            rs_vec = collect(rs)
            thetas_vec = collect(thetas)
            op_pair_arr = op_func(factor * qpairbasis, rs_vec, thetas_vec)
            op_block_arr = op_func(factor * qblockbasis, rs_vec, thetas_vec)
            @test op_pair_arr isa GaussianUnitary && op_block_arr isa GaussianUnitary
            @test op_func(Array, factor * qpairbasis, rs_vec, thetas_vec) isa GaussianUnitary
            @test op_func(SVector{factor * 2 * nmodes}, SMatrix{factor * 2 * nmodes, factor * 2 * nmodes}, factor * qpairbasis, rs_vec, thetas_vec) isa GaussianUnitary
            @test all.(isequal(op_func(factor * qblockbasis, rs_vec, thetas_vec).disp, changebasis(QuadBlockBasis, op_pair_arr).disp))
            @test all.(isequal(op_func(factor * qblockbasis, rs_vec, thetas_vec).symplectic, changebasis(QuadBlockBasis, op_pair_arr).symplectic))
        end
    end

    for (name, op_func, factor) in [("phase-shift", phaseshift, 1), ("beamsplitter", beamsplitter, 2)]
        @testset "Symbolic $name operator" begin
            @variables theta
            op = op_func(factor * qpairbasis, theta)
            @test op isa GaussianUnitary
            @test op_func(Array, factor * qpairbasis, theta) isa GaussianUnitary
            @test op_func(SVector{factor * 2 * nmodes}, SMatrix{factor * 2 * nmodes, factor * 2 * nmodes}, factor * qpairbasis, theta) isa GaussianUnitary
            @test all.(isequal(op_func(factor * qblockbasis, theta).disp, changebasis(QuadBlockBasis, op).disp))
            @test all.(isequal(op_func(factor * qblockbasis, theta).symplectic, changebasis(QuadBlockBasis, op).symplectic))
            @variables thetas[1:nmodes]
            thetas_vec = collect(thetas)
            op_arr = op_func(factor * qpairbasis, thetas_vec)
            op_block_arr = op_func(factor * qblockbasis, thetas_vec)
            @test op_arr isa GaussianUnitary && op_block_arr isa GaussianUnitary
            @test op_func(Array, factor * qpairbasis, thetas_vec) isa GaussianUnitary
            @test op_func(SVector{factor * 2 * nmodes}, SMatrix{factor * 2 * nmodes, factor * 2 * nmodes}, factor * qpairbasis, thetas_vec) isa GaussianUnitary
            @test all.(isequal(op_func(factor * qblockbasis, thetas_vec).disp, changebasis(QuadBlockBasis, op_arr).disp))
            @test all.(isequal(op_func(factor * qblockbasis, thetas_vec).symplectic, changebasis(QuadBlockBasis, op_arr).symplectic))
        end
    end
end

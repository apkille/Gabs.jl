@testitem "Unitaries" begin
    using Gabs
    using StaticArrays

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)

    @testset "displacement operator" begin
        alpha = rand(ComplexF64)
        alphas = rand(ComplexF64, nmodes)
        op_pair = displace(qpairbasis, alpha)
        op_block = displace(qblockbasis, alpha)
        op, op_array, op_static = displace(qpairbasis, alpha), displace(Array, qpairbasis, alpha), displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, alpha)
        op_static_array, op_static_vec = displace(SArray, qpairbasis,alpha), displace(SVector, SMatrix, qpairbasis, alpha)
        @test op_pair isa GaussianUnitary && op_block isa GaussianUnitary
        @test op isa GaussianUnitary && op_array isa GaussianUnitary && op_static isa GaussianUnitary
        @test op_static_array isa GaussianUnitary && op_static_vec isa GaussianUnitary
        @test displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, alpha) isa GaussianUnitary
        @test op_pair == changebasis(QuadPairBasis, op_block) && op_block == changebasis(QuadBlockBasis, op_pair)
        @test op_pair == changebasis(QuadPairBasis, op_pair) && op_block == changebasis(QuadBlockBasis, op_block)
        @test displace(qblockbasis, alpha) == changebasis(QuadBlockBasis, op_pair)
        @test displace(qblockbasis, alphas) == changebasis(QuadBlockBasis, displace(qpairbasis, alphas))
        @test issymplectic(qpairbasis, op_pair.symplectic, atol = 1e-4)
        @test isgaussian(op_pair, atol = 1e-4)
        @test op_pair.ħ == 2 && op_block.ħ == 2
    end

    @testset "squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        rs, thetas = rand(Float64, nmodes), rand(Float64, nmodes)
        op_pair = squeeze(qpairbasis, r, theta)
        op_block = squeeze(qblockbasis, r, theta)
        op, op_array, op_static = squeeze(2*qpairbasis, r, theta), squeeze(Array, 2*qpairbasis, r, theta), squeeze(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, r, theta)
        op_static_array, op_static_vec = squeeze(SArray, 2*qpairbasis, r, theta), squeeze(SVector, SMatrix, 2*qpairbasis, r, theta)
        @test op_pair isa GaussianUnitary && op_block isa GaussianUnitary
        @test op isa GaussianUnitary && op_array isa GaussianUnitary && op_static isa GaussianUnitary
        @test op_static_array isa GaussianUnitary && op_static_vec isa GaussianUnitary
        @test op_pair == changebasis(QuadPairBasis, op_block) && op_block == changebasis(QuadBlockBasis, op_pair)
        @test op_pair == changebasis(QuadPairBasis, op_pair) && op_block == changebasis(QuadBlockBasis, op_block)
        @test squeeze(qblockbasis, r, theta) == changebasis(QuadBlockBasis, op_pair)
        @test squeeze(qblockbasis, rs, thetas) == changebasis(QuadBlockBasis, squeeze(qpairbasis, rs, thetas))
        @test issymplectic(qpairbasis, op_pair.symplectic, atol = 1e-4)
        @test isgaussian(op_pair, atol = 1e-4)
        @test op_pair.ħ == 2 && op_block.ħ == 2
    end

    @testset "two-mode squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        rs, thetas = rand(Float64, nmodes), rand(Float64, nmodes)
        op, op_array, op_static = twosqueeze(2*qpairbasis, r, theta), twosqueeze(Array, 2*qpairbasis, r, theta), twosqueeze(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, r, theta)
        op_static_array, op_static_vec = twosqueeze(SArray, 2*qpairbasis, r, theta), twosqueeze(SVector, SMatrix, 2*qpairbasis, r, theta)
        @test op isa GaussianUnitary && op_array isa GaussianUnitary && op_static isa GaussianUnitary
        @test op_static_array isa GaussianUnitary && op_static_vec isa GaussianUnitary
        @test twosqueeze(2*qblockbasis, r, theta) == changebasis(QuadBlockBasis, op)
        @test twosqueeze(2*qblockbasis, rs, thetas) == changebasis(QuadBlockBasis, twosqueeze(2*qpairbasis, rs, thetas))
        @test op.ħ == 2 && op_array.ħ == 2 && op_static.ħ == 2
    end

    @testset "phase-shift operator" begin
        theta = rand(Float64)
        thetas = rand(Float64, nmodes)
        op, op_array, op_static = phaseshift(qpairbasis, theta), phaseshift(Array, qpairbasis, theta), phaseshift(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, theta)
        op_static_array, op_static_vec = phaseshift(SArray, 2*qpairbasis, theta), phaseshift(SVector, SMatrix, 2*qpairbasis, theta)
        @test op isa GaussianUnitary && op_array isa GaussianUnitary && op_static isa GaussianUnitary
        @test op_static_array isa GaussianUnitary && op_static_vec isa GaussianUnitary
        @test phaseshift(qblockbasis, theta) == changebasis(QuadBlockBasis, op)
        @test phaseshift(qblockbasis, thetas) == changebasis(QuadBlockBasis, phaseshift(qpairbasis, thetas))
        @test op.ħ == 2 && op_array.ħ == 2 && op_static.ħ == 2
    end

    @testset "beamsplitter operator" begin
        theta = rand(Float64)
        thetas = rand(Float64, nmodes)
        op, op_array, op_static = beamsplitter(2*qpairbasis, theta), beamsplitter(Array, 2*qpairbasis, theta), beamsplitter(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, theta)
        op_static_array, op_static_vec = beamsplitter(SArray, 2*qpairbasis, theta), beamsplitter(SVector, SMatrix, 2*qpairbasis, theta)
        @test op isa GaussianUnitary && op_array isa GaussianUnitary && op_static isa GaussianUnitary
        @test op_static_array isa GaussianUnitary && op_static_vec isa GaussianUnitary
        @test beamsplitter(2*qblockbasis, theta) == changebasis(QuadBlockBasis, op)
        @test beamsplitter(2*qblockbasis, thetas) == changebasis(QuadBlockBasis, beamsplitter(2*qpairbasis, thetas))
        @test op.ħ == 2 && op_array.ħ == 2 && op_static.ħ == 2
    end

    @testset "tensor products" begin
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(qpairbasis, alpha1), displace(qpairbasis, alpha2)
        ds = tensor(d1, d2)
        @test ds isa GaussianUnitary
        @test ds == d1 ⊗ d2
        @test isapprox(ds, d1 ⊗ d2, atol = 1e-10)
        @test tensor(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, d1, d2) isa GaussianUnitary

        r, theta = rand(Float64), rand(Float64)
        p = phaseshift(qpairbasis, theta)
        @test tensor(p, tensor(d1, d2)) == p ⊗ d1 ⊗ d2

        p_block = phaseshift(qblockbasis, theta)
        p_blocks = phaseshift(2*qblockbasis, repeat([theta], 2*nmodes))
        @test p_block ⊗ p_block == p_blocks

        dstatic, dstatic1 = displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, alpha1), displace(SVector, SMatrix, qpairbasis, alpha1)
        tpstatic, tpstatic1 = dstatic ⊗ dstatic ⊗ dstatic,  dstatic1 ⊗ dstatic1 ⊗ dstatic1
        @test tpstatic.disp isa SVector{6*nmodes} && tpstatic1.disp isa SVector{6*nmodes}
        @test tpstatic.symplectic isa SMatrix{6*nmodes,6*nmodes} && tpstatic1.symplectic isa SMatrix{6*nmodes,6*nmodes}
        tp, tp1 = dstatic ⊗ d1 ⊗ dstatic, dstatic1 ⊗ d1 ⊗ dstatic1
        @test tp.disp isa SVector && tp1.disp isa SVector
        @test tp.symplectic isa SMatrix && tp1.symplectic isa SMatrix
    end

    @testset "actions" begin
        alpha = rand(ComplexF64)
        d = displace(qpairbasis, alpha)
        v = vacuumstate(qpairbasis)
        c = coherentstate(qpairbasis, alpha)
        @test d * v == c
        @test apply!(v, d) == c
        
        v1, v2 = vacuumstate(qpairbasis), vacuumstate(qpairbasis)
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(qpairbasis, alpha1), displace(qpairbasis, alpha2)
        c1, c2 = coherentstate(qpairbasis, alpha1), coherentstate(qpairbasis, alpha2)
        @test (d1 ⊗ d2) * (v1 ⊗ v2) == c1 ⊗ c2
    end
end
@testitem "Unitaries" begin
    import Gabs: _changebasis
    using Gabs
    using StaticArrays
    
    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)

    @testset "displacement operator" begin
        alpha = rand(ComplexF64)
        alphas = rand(ComplexF64, nmodes)
        op = displace(qpairbasis, alpha)
        @test op isa GaussianUnitary
        @test displace(Array, qpairbasis, alpha) isa GaussianUnitary
        @test displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, alpha) isa GaussianUnitary
        @test displace(qblockbasis, alpha) == _changebasis(op, QuadBlockBasis)
        @test displace(qblockbasis, alphas) == _changebasis(displace(qpairbasis, alphas), QuadBlockBasis)
    end

    @testset "squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        rs, thetas = rand(Float64, nmodes), rand(Float64, nmodes)
        op = squeeze(qpairbasis, r, theta)
        @test op isa GaussianUnitary
        @test squeeze(Array, qpairbasis, r, theta) isa GaussianUnitary
        @test squeeze(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, r, theta) isa GaussianUnitary
        @test squeeze(qblockbasis, r, theta) == _changebasis(op, QuadBlockBasis)
        @test squeeze(qblockbasis, rs, thetas) == _changebasis(squeeze(qpairbasis, rs, thetas), QuadBlockBasis)
    end

    @testset "two-mode squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        rs, thetas = rand(Float64, nmodes), rand(Float64, nmodes)
        op = twosqueeze(2*qpairbasis, r, theta)
        @test op isa GaussianUnitary
        @test twosqueeze(Array, 2*qpairbasis, r, theta) isa GaussianUnitary
        @test twosqueeze(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, r, theta) isa GaussianUnitary
        @test twosqueeze(2*qblockbasis, r, theta) == _changebasis(op, QuadBlockBasis)
        @test twosqueeze(2*qblockbasis, rs, thetas) == _changebasis(twosqueeze(2*qpairbasis, rs, thetas), QuadBlockBasis)
    end

    @testset "phase-shift operator" begin
        theta = rand(Float64)
        thetas = rand(Float64, nmodes)
        op = phaseshift(qpairbasis, theta)
        @test op isa GaussianUnitary
        @test phaseshift(Array, qpairbasis, theta) isa GaussianUnitary
        @test phaseshift(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, theta) isa GaussianUnitary
        @test phaseshift(qblockbasis, theta) == _changebasis(op, QuadBlockBasis)
        @test phaseshift(qblockbasis, thetas) == _changebasis(phaseshift(qpairbasis, thetas), QuadBlockBasis)
    end

    @testset "beamsplitter operator" begin
        theta = rand(Float64)
        thetas = rand(Float64, nmodes)
        op = beamsplitter(2*qpairbasis, theta)
        @test op isa GaussianUnitary
        @test beamsplitter(Array, 2*qpairbasis, theta) isa GaussianUnitary
        @test beamsplitter(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, theta) isa GaussianUnitary
        @test beamsplitter(2*qblockbasis, theta) == _changebasis(op, QuadBlockBasis)
        @test beamsplitter(2*qblockbasis, thetas) == _changebasis(beamsplitter(2*qpairbasis, thetas), QuadBlockBasis)
    end

    @testset "tensor products" begin
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(qpairbasis, alpha1), displace(qpairbasis, alpha2)
        ds = tensor(d1, d2)
        @test ds isa GaussianUnitary
        @test ds == d1 ⊗ d2
        @test isapprox(ds, d1 ⊗ d2)
        @test tensor(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, d1, d2) isa GaussianUnitary

        r, theta = rand(Float64), rand(Float64)
        p = phaseshift(qpairbasis, theta)
        @test tensor(p, tensor(d1, d2)) == p ⊗ d1 ⊗ d2

        p_block = phaseshift(qblockbasis, theta)
        p_blocks = phaseshift(2*qblockbasis, repeat([theta], 2*nmodes))
        @test p_block ⊗ p_block == p_blocks

        dstatic = displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, alpha1)
        tpstatic = dstatic ⊗ dstatic ⊗ dstatic
        @test tpstatic.disp isa SVector{6*nmodes}
        @test tpstatic.symplectic isa SMatrix{6*nmodes,6*nmodes}
        tp = dstatic ⊗ d1 ⊗ dstatic
        @test tp.disp isa Vector
        @test tp.symplectic isa Matrix
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
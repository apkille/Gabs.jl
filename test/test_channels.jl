@testitem "Channels" begin
    import Gabs: _changebasis
    using Gabs
    using StaticArrays

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)
    noise = rand(2*nmodes, 2*nmodes)
    noise_ds = [noise zeros(2*nmodes,2*nmodes); zeros(2*nmodes,2*nmodes) noise]
    T = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(2*nmodes), j in Base.OneTo(2*nmodes)
        if (j == 2*i-1) || (j + 2*nmodes == 2*i)
            T[i,j] = 1
        end
    end
    T_ds = zeros(4*nmodes, 4*nmodes)
    @inbounds for i in Base.OneTo(4*nmodes), j in Base.OneTo(4*nmodes)
        if (j == 2*i-1) || (j + 2*(2*nmodes) == 2*i)
            T_ds[i,j] = 1
        end
    end

    @testset "displacement operator" begin
        alpha = rand(ComplexF64)
        op = displace(qpairbasis, alpha, noise)
        @test op isa GaussianChannel
        @test displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, alpha, noise) isa GaussianChannel
        @test displace(Array, qpairbasis, alpha, noise) isa GaussianChannel
        @test displace(qblockbasis, alpha, T*noise*transpose(T)) == _changebasis(op, QuadBlockBasis)
    end

    @testset "squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        op = squeeze(qpairbasis, r, theta, noise)
        @test op isa GaussianChannel
        @test squeeze(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, r, theta, noise) isa GaussianChannel
        @test squeeze(Array, qpairbasis, r, theta, noise) isa GaussianChannel
        @test squeeze(qblockbasis, r, theta, T*noise*transpose(T)) == _changebasis(op, QuadBlockBasis)
    end

    @testset "two-mode squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        op = twosqueeze(2*qpairbasis, r, theta, noise_ds)
        @test op isa GaussianChannel
        @test twosqueeze(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, r, theta, noise_ds) isa GaussianChannel
        @test twosqueeze(Array, 2*qpairbasis, r, theta, noise_ds) isa GaussianChannel
        @test twosqueeze(2*qblockbasis, r, theta, T_ds*noise_ds*transpose(T_ds)) == _changebasis(op, QuadBlockBasis)
    end

    @testset "phase-shift operator" begin
        theta = rand(Float64)
        op = phaseshift(qpairbasis, theta, noise)
        @test op isa GaussianChannel
        @test phaseshift(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, theta, noise) isa GaussianChannel
        @test phaseshift(Array, qpairbasis, theta, noise) isa GaussianChannel
        @test phaseshift(qblockbasis, theta, T*noise*transpose(T)) == _changebasis(op, QuadBlockBasis)
    end

    @testset "beamsplitter operator" begin
        theta = rand(Float64)
        op = beamsplitter(2*qpairbasis, theta, noise_ds)
        @test op isa GaussianChannel
        @test beamsplitter(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, theta, noise_ds) isa GaussianChannel
        @test beamsplitter(Array, 2*qpairbasis, theta, noise_ds) isa GaussianChannel
        @test beamsplitter(2*qblockbasis, theta, T_ds*noise_ds*transpose(T_ds)) == _changebasis(op, QuadBlockBasis)
    end

    @testset "attenuator channel" begin
        theta = rand(Float64)
        n = rand(Int64)
        op = attenuator(qpairbasis, theta, n)
        @test op isa GaussianChannel
        @test attenuator(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, theta, n) isa GaussianChannel
        @test attenuator(Array, qpairbasis, theta, n) isa GaussianChannel
        @test attenuator(qblockbasis, theta, n) == _changebasis(op, QuadBlockBasis)
    end

    @testset "amplifier channel" begin
        r = rand(Float64)
        n = rand(Int64)
        op = amplifier(qpairbasis, r, n)
        @test op isa GaussianChannel
        @test amplifier(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, r, n) isa GaussianChannel
        @test amplifier(Array, qpairbasis, r, n) isa GaussianChannel
        @test amplifier(qblockbasis, r, n) == _changebasis(op, QuadBlockBasis)
    end
    
    @testset "tensor products" begin
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(qpairbasis, alpha1, noise), displace(qpairbasis, alpha2, noise)
        ds = tensor(d1, d2)
        @test ds isa GaussianChannel
        @test ds == d1 ⊗ d2
        @test tensor(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, d1, d2) isa GaussianChannel

        r, theta = rand(Float64), rand(Float64)
        p = phaseshift(qpairbasis, theta, noise)
        @test tensor(tensor(p, d1), d2) == p ⊗ d1 ⊗ d2


        dstatic = displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, alpha1, noise)
        tpstatic = dstatic ⊗ dstatic ⊗ dstatic
        @test tpstatic.disp isa SVector{6*nmodes}
        @test tpstatic.transform isa SMatrix{6*nmodes,6*nmodes}
        @test tpstatic.noise isa SMatrix{6*nmodes,6*nmodes}
        tp = dstatic ⊗ d1 ⊗ dstatic
        @test tp.disp isa Vector
        @test tp.transform isa Matrix
        @test tp.noise isa Matrix
    end

    @testset "actions" begin
        z = zeros(2*nmodes,2*nmodes)
        alpha = rand(ComplexF64)
        d = displace(qpairbasis, alpha, z)
        v = vacuumstate(qpairbasis)
        c = coherentstate(qpairbasis, alpha)
        @test d * v == c
        @test isapprox(d * v, c)
        @test apply!(v, d) == c

        v1, v2 = vacuumstate(qpairbasis), vacuumstate(qpairbasis)
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(qpairbasis, alpha1, z), displace(qpairbasis, alpha2, z)
        c1, c2 = coherentstate(qpairbasis, alpha1), coherentstate(qpairbasis, alpha2)
        @test (d1 ⊗ d2) * (v1 ⊗ v2) == c1 ⊗ c2
    end
end
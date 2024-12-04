@testitem "Quadrature pair basis" begin
    using Gabs
    using StaticArrays

    basis1 = QuadPairBasis(1)
    basis2 = QuadPairBasis(2)

    noise2 = rand(Float64, (2,2))
    noise2_ds = [noise2 zeros(2,2); zeros(2,2) noise2]
    noise4 = rand(Float64, (4,4))

    @testset "displacement operator" begin
        alpha = rand(ComplexF64)
        @test displace(basis1, alpha, noise2) isa GaussianChannel
        @test displace(SVector{2}, SMatrix{2,2}, basis1, alpha, noise2) isa GaussianChannel
        @test displace(Array, basis1, alpha, noise2) isa GaussianChannel
    end

    @testset "squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        @test squeeze(basis1, r, theta, noise2) isa GaussianChannel
        @test squeeze(SVector{2}, SMatrix{2,2}, basis1, r, theta, noise2) isa GaussianChannel
        @test squeeze(Array, basis1, r, theta, noise2) isa GaussianChannel
    end

    @testset "two-mode squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        @test twosqueeze(basis2, r, theta, noise4) isa GaussianChannel
        @test twosqueeze(SVector{4}, SMatrix{4,4}, basis2, r, theta, noise4) isa GaussianChannel
        @test twosqueeze(Array, basis2, r, theta, noise4) isa GaussianChannel
    end

    @testset "phase-shift operator" begin
        theta = rand(Float64)
        @test phaseshift(basis1, theta, noise2) isa GaussianChannel
        @test phaseshift(SVector{2}, SMatrix{2,2}, basis1, theta, noise2) isa GaussianChannel
        @test phaseshift(Array, basis1, theta, noise2) isa GaussianChannel
    end

    @testset "beamsplitter operator" begin
        theta = rand(Float64)
        @test beamsplitter(basis2, theta, noise4) isa GaussianChannel
        @test beamsplitter(SVector{4}, SMatrix{4,4}, basis2, theta, noise4) isa GaussianChannel
        @test beamsplitter(Array, basis2, theta, noise4) isa GaussianChannel
    end

    @testset "attenuator channel" begin
        theta = rand(Float64)
        n = rand(Int64)
        @test attenuator(basis1, theta, n) isa GaussianChannel
        @test attenuator(SVector{2}, SMatrix{2,2}, basis1, theta, n) isa GaussianChannel
        @test attenuator(Array, basis1, theta, n) isa GaussianChannel
    end

    @testset "amplifier channel" begin
        r = rand(Float64)
        n = rand(Int64)
        @test amplifier(basis1, r, n) isa GaussianChannel
        @test amplifier(SVector{2}, SMatrix{2,2}, basis1, r, n) isa GaussianChannel
        @test amplifier(Array, basis1, r, n) isa GaussianChannel
    end
    
    @testset "tensor products" begin
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(basis1, alpha1, noise2), displace(basis1, alpha2, noise2)
        ds = tensor(d1, d2)
        @test ds isa GaussianChannel
        @test ds == d1 ⊗ d2
        @test tensor(SVector{4}, SMatrix{4,4}, d1, d2) isa GaussianChannel

        r, theta = rand(Float64), rand(Float64)
        p = phaseshift(basis1, theta, noise2)
        @test tensor(tensor(p, d1), d2) == p ⊗ d1 ⊗ d2


        dstatic = displace(SVector{2}, SMatrix{2,2}, basis1, alpha1, noise2)
        tpstatic = dstatic ⊗ dstatic ⊗ dstatic
        @test tpstatic.disp isa SVector{6}
        @test tpstatic.transform isa SMatrix{6,6}
        @test tpstatic.noise isa SMatrix{6,6}
        tp = dstatic ⊗ d1 ⊗ dstatic
        @test tp.disp isa Vector
        @test tp.transform isa Matrix
        @test tp.noise isa Matrix
    end

    @testset "actions" begin
        z = zeros(2,2)
        alpha = rand(ComplexF64)
        d = displace(basis1, alpha, z)
        v = vacuumstate(basis1)
        c = coherentstate(basis1, alpha)
        @test d * v == c
        @test isapprox(d * v, c)
        @test apply!(v, d) == c

        v1, v2 = vacuumstate(basis1), vacuumstate(basis1)
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(basis1, alpha1, z), displace(basis1, alpha2, z)
        c1, c2 = coherentstate(basis1, alpha1), coherentstate(basis1, alpha2)
        @test (d1 ⊗ d2) * (v1 ⊗ v2) == c1 ⊗ c2
    end
end
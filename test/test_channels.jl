@testitem "Channels" begin
    using Gabs
    using StaticArrays

    noise2 = rand(Float64, (2,2))
    noise2_ds = [noise2 zeros(2,2); zeros(2,2) noise2]
    noise4 = rand(Float64, (4,4))

    @testset "displacement operator" begin
        alpha = rand(ComplexF64)
        @test displace(alpha, noise2) isa GaussianChannel
        @test displace(SVector{2}, SMatrix{2,2}, alpha, noise2) isa GaussianChannel
    end

    @testset "squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        @test squeeze(r, theta, noise2) isa GaussianChannel
        @test squeeze(SVector{2}, SMatrix{2,2}, r, theta, noise2) isa GaussianChannel
    end

    @testset "two-mode squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        @test twosqueeze(r, theta, noise4) isa GaussianChannel
        @test twosqueeze(SVector{4}, SMatrix{4,4}, r, theta, noise4) isa GaussianChannel
    end

    @testset "phase-shift operator" begin
        theta = rand(Float64)
        @test phaseshift(theta, noise2) isa GaussianChannel
        @test phaseshift(SVector{2}, SMatrix{2,2}, theta, noise2) isa GaussianChannel
    end

    @testset "beamsplitter operator" begin
        theta = rand(Float64)
        @test beamsplitter(theta, noise4) isa GaussianChannel
        @test beamsplitter(SVector{4}, SMatrix{4,4}, theta, noise4) isa GaussianChannel
    end

    @testset "attenuator channel" begin
        theta = rand(Float64)
        n = rand(Int64)
        @test attenuator(theta, n) isa GaussianChannel
        @test attenuator(SVector{2}, SMatrix{2,2}, theta, n) isa GaussianChannel
    end

    @testset "attenuator channel" begin
        theta = rand(Float64)
        n = rand(Int64)
        @test attenuator(theta, n) isa GaussianChannel
        @test attenuator(SVector{2}, SMatrix{2,2}, theta, n) isa GaussianChannel
    end
    
    @testset "tensor products" begin
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(alpha1, noise2), displace(alpha2, noise2)
        ds = tensor(d1, d2)
        @test ds isa GaussianChannel
        @test ds == d1 ⊗ d2
        @test tensor(SVector{4}, SMatrix{4,4}, d1, d2) isa GaussianChannel

        r, theta = rand(Float64), rand(Float64)
        p = phaseshift(theta, noise2)
        @test tensor(tensor(p, d1), d2) == p ⊗ d1 ⊗ d2

        s1, s2 = squeeze(r, theta, noise2), squeeze(r, theta, noise2)
        @testset s1 ⊗ s2 == twosqueeze(r, theta, noise2_ds)
    end

    @testset "actions" begin
        z = zeros(2,2)
        alpha = rand(ComplexF64)
        d = displace(alpha, z)
        v = vacuumstate()
        c = coherentstate(alpha)
        @test apply(v, d) == c
        @test apply!(v, d) == c

        v1, v2 = vacuumstate(), vacuumstate()
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(alpha1, z), displace(alpha2, z)
        c1, c2 = coherentstate(alpha1), coherentstate(alpha2)
        @test apply(v1 ⊗ v2, d1 ⊗ d2) == c1 ⊗ c2
    end
end
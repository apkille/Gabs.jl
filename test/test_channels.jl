@testitem "Channels" begin
    using Gabs, StaticArrays

    @testset "displacement operator" begin
        @test_nowarn displace(1.0+pi*im)
        @test_nowarn displace(SVector{2}, SMatrix{2,2}, 1.0+pi*im)
    end

    @testset "squeeze operator" begin
        @test_nowarn squeeze(1.0, 2.0)
        @test_nowarn squeeze(SVector{2}, SMatrix{2,2}, 1.0, 2.0)
    end

    @testset "two-mode squeeze operator" begin
        @test_nowarn twosqueeze(1.0, 2.0)
        @test_nowarn twosqueeze(SVector{4}, SMatrix{4,4}, 1.0, 2.0)
    end

    @testset "phase-shift operator" begin
        @test_nowarn phaseshift(2.0)
        @test_nowarn phaseshift(SVector{2}, SMatrix{2,2}, 2.0)
    end

    @testset "beamsplitter operator" begin
        @test_nowarn beamsplitter(0.25)
        @test_nowarn beamsplitter(SVector{4}, SMatrix{4,4}, 0.25)
    end

    @testset "direct sums" begin
        d1, d2 = displace(1.0+pi*im), displace(4.0-3.0*im)
        @test directsum(d1, d2) == d1 ⊕ d2
        @test_nowarn directsum(SVector{4}, SMatrix{4,4}, d1, d2)
        p = phaseshift(-9.0)
        @test directsum(p, directsum(d1, d2)) == p ⊕ d1 ⊕ d2
        s1, s2 = squeeze(1.0, 2.0), squeeze(1.0, 2.0)
        @testset s1 ⊕ s2 == twosqueeze(1.0, 2.0)
    end

    @testset "actions" begin
        d = displace(1.0-im)
        v = vacuumstate()
        c = coherentstate(1.0-im)
        @test apply(v, d) == c
        @test apply!(v, d) == c
        ts = twosqueeze(1.0, 2.0)
        v1, v2 = vacuumstate(), vacuumstate()
        d1, d2 = displace(1.0-im), displace(3.0+9.0im)
        c1, c2 = coherentstate(1.0-im), coherentstate(3.0+9.0im)
        @test apply(v1 ⊕ v2, d1 ⊕ d2) == c1 ⊕ c2
    end
end
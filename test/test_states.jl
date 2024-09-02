@testitem "States" begin
    using Gabs, StaticArrays

    @testset "vacuum states" begin
        @test_nowarn vacuumstate()
        @test_nowarn vacuumstate(SVector{2}, SMatrix{2,2})
    end

    @testset "thermal states" begin
        @test_nowarn thermalstate(2)
        @test_nowarn thermalstate(SVector{2}, SMatrix{2,2}, 2)
    end

    @testset "coherent states" begin
        @test_nowarn coherentstate(1.0+pi*im)
        @test_nowarn coherentstate(SVector{2}, SMatrix{2,2}, 1.0+pi*im)
    end

    @testset "squeezed states" begin
        @test_nowarn squeezedstate(1.0, 2.0)
        @test_nowarn squeezedstate(SVector{2}, SMatrix{2,2}, 1.0, 2.0)
    end

    @testset "epr states" begin
        @test_nowarn eprstate(1.0, 2.0)
        @test_nowarn eprstate(SVector{4}, SMatrix{4,4}, 1.0, 2.0)
    end

    @testset "direct sums" begin
        v1, v2 = vacuumstate(), vacuumstate()
        @test directsum(v1, v2) == v1 ⊕ v2
        @test_nowarn directsum(SVector{4}, SMatrix{4,4}, v1, v2)
        c = coherentstate(1.0 + pi*im)
        @test directsum(c, directsum(v1, v2)) == c ⊕ v1 ⊕ v2
        s1, s2 = squeezedstate(1.0, 2.0), squeezedstate(1.0, 2.0)
        @testset s1 ⊕ s2 == eprstate(1.0, 2.0)
    end
end
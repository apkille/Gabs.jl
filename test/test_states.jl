@testitem "States" begin
    using Gabs, StaticArrays

    @testset "vacuum states" begin
        @test vacuumstate() isa GaussianState
        @test vacuumstate(SVector{2}, SMatrix{2,2}) isa GaussianState
    end

    @testset "thermal states" begin
        n = rand(Int64)
        @test thermalstate(n) isa GaussianState
        @test thermalstate(SVector{2}, SMatrix{2,2}, n) isa GaussianState
    end

    @testset "coherent states" begin
        alpha = rand(ComplexF64)
        @test coherentstate(alpha) isa GaussianState
        @test coherentstate(SVector{2}, SMatrix{2,2}, alpha) isa GaussianState
    end

    @testset "squeezed states" begin
        r, theta = rand(Float64), rand(Float64)
        @test squeezedstate(r, theta) isa GaussianState
        @test squeezedstate(SVector{2}, SMatrix{2,2}, r, theta) isa GaussianState
    end

    @testset "epr states" begin
        r, theta = rand(Float64), rand(Float64)
        @test eprstate(r, theta) isa GaussianState
        @test eprstate(SVector{4}, SMatrix{4,4}, r, theta) isa GaussianState
    end

    @testset "direct sums" begin
        v1, v2 = vacuumstate(), vacuumstate()
        ds = directsum(v1, v2)
        @test ds isa GaussianState
        @test directsum(SVector{4}, SMatrix{4,4}, v1, v2) isa GaussianState
        @test ds == v1 ⊕ v2

        alpha = rand(ComplexF64)
        c = coherentstate(alpha)
        @test directsum(c, directsum(v1, v2)) == c ⊕ v1 ⊕ v2

        r, theta = rand(Float64), rand(Float64)
        s1, s2 = squeezedstate(r, theta), squeezedstate(r, theta)
        @testset s1 ⊕ s2 == eprstate(r, theta)
    end
end
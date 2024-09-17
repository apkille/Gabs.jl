@testitem "Throws" begin
    using Gabs
    using CairoMakie

    @testset "type throws" begin
        @test_throws DimensionMismatch GaussianState([1.0, 2.0, 3.0], [3.0 4.0; 5.0 6.0])
        @test_throws DimensionMismatch GaussianChannel([1.0, 2.0, 3.0], [3.0 4.0; 5.0 6.0], [3.0 4.0; 5.0 6.0])
        @test_throws DimensionMismatch GaussianChannel([1.0, 2.0], [3.0 4.0; 5.0 6.0], [3.0 4.0 4.0; 5.0 6.0 4.0])
    end

    @testset "action throws" begin
        v = vacuumstate()
        ts = twosqueeze(rand(), rand())
        @test_throws DimensionMismatch ts * v
        @test_throws DimensionMismatch apply!(v, ts)
    end

    @testset "plot extension throws" begin
        ts = twosqueeze(rand(), rand())
        @test_throws ArgumentError Makie.heatmap(collect(-3.0:0.25:3.0), collect(-3.0:0.25:3.0), ts)
    end
end
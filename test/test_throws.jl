@testitem "Throws" begin
    using Gabs
    using CairoMakie

    basis1 = QuadPairBasis(1)
    @testset "type throws" begin
        @test_throws DimensionMismatch GaussianState(basis1, [1.0, 2.0, 3.0], [3.0 4.0; 5.0 6.0])
        @test_throws DimensionMismatch GaussianChannel(basis1, [1.0, 2.0, 3.0], [3.0 4.0; 5.0 6.0], [3.0 4.0; 5.0 6.0])
        @test_throws DimensionMismatch GaussianChannel(basis1, [1.0, 2.0], [3.0 4.0; 5.0 6.0], [3.0 4.0 4.0; 5.0 6.0 4.0])
        vac = vacuumstate(basis1)
        vacs = vac ⊗ vac ⊗ vac ⊗ vac
        @test_throws DimensionMismatch Generaldyne(vacs, vac ⊗ vac, [2, 4, 5])
    end

    @testset "action throws" begin
        v = vacuumstate(basis1)
        ts = twosqueeze(2*basis1, rand(), rand())
        @test_throws DimensionMismatch ts * v
        @test_throws DimensionMismatch apply!(v, ts)
    end

    @testset "plot extension throws" begin
        ts = twosqueeze(2*basis1, rand(), rand())
        @test_throws ArgumentError Makie.heatmap(collect(-3.0:0.25:3.0), collect(-3.0:0.25:3.0), ts)
    end
end
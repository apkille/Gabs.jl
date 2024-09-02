@testitem "Throws" begin
    using Gabs

    @testset "type throws" begin
        @test_throws DimensionMismatch GaussianState([1.0, 2.0, 3.0], [3.0 4.0; 5.0 6.0])
        @test_throws DimensionMismatch GaussianChannel([1.0, 2.0, 3.0], [3.0 4.0; 5.0 6.0], [3.0 4.0; 5.0 6.0])
        @test_throws DimensionMismatch GaussianChannel([1.0, 2.0], [3.0 4.0; 5.0 6.0], [3.0 4.0 4.0; 5.0 6.0 4.0])
    end

    @testset "action throws" begin
        v = vacuumstate()
        ts = twosqueeze(1.0, 4.0)
        @test_throws DimensionMismatch apply(v, ts)
        @test_throws DimensionMismatch apply!(v, ts)
    end
end
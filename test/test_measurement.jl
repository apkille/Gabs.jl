@testitem "Measurements" begin
    using Gabs

    @testset "displacement operator" begin
        alpha = rand(ComplexF64)
        @test displace(alpha, noise2) isa GaussianChannel
        @test displace(SVector{2}, SMatrix{2,2}, alpha, noise2) isa GaussianChannel
    end
end
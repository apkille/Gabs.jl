@testitem "Measurements" begin
    using Gabs
    using StaticArrays

    @testset "random types" begin
        nmodes = rand(1:20)
        rs = randstate(nmodes)
        rc = randchannel(nmodes)
        @test rc isa GaussianChannel
        @test rc * rs isa GaussianState

        rsfloat32 = randstate(AbstractArray{Float32}, nmodes)
        rcfloat32 = randchannel(AbstractArray{Float32}, nmodes)
        @test rcfloat32 isa GaussianChannel
        @test rcfloat32 * rsfloat32 isa GaussianState

        rsstatic = randstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, nmodes)
        rcstatic = randchannel(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, nmodes)
        @test rcstatic isa GaussianChannel
        @test rcstatic * rsstatic isa GaussianState
    end
end
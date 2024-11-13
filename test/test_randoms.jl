@testitem "Measurements" begin
    using Gabs
    using StaticArrays
    using LinearAlgebra: eigvals, adjoint

    @testset "random utils" begin
        nmodes = rand(1:20)
        U = Gabs._rand_unitary(nmodes)
        @test isapprox(adjoint(U), inv(U))

        O = Gabs._rand_orthogonal_symplectic(nmodes)
        @test isapprox(O', inv(O))
        Omega = symplecticform(nmodes)
        @test isapprox(O' * Omega * O, Omega)
    end

    @testset "random types" begin
        nmodes = rand(1:20)
        rs = randstate(nmodes)
        rc = randchannel(nmodes)
        @test rc isa GaussianChannel
        @test rc * rs isa GaussianState
        Omega = symplecticform(nmodes)
        @test all(i -> (i >= 0) || isapprox(i, 0.0, atol = 1e-5), real(eigvals(rs.covar .+ im*Omega)))

        rspure = randstate(nmodes, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rspure.covar .+ im*Omega)))
        @test isapprox(purity(rspure), 1.0, atol = 1e-5)

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
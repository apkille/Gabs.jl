@testitem "Measurements" begin
    using Gabs
    using StaticArrays
    using LinearAlgebra: eigvals, adjoint

    @testset "random utils" begin
        nmodes = rand(1:20)
        U = Gabs._rand_unitary(nmodes)
        @test isapprox(adjoint(U), inv(U), atol = 1e-5)

        O = Gabs._rand_orthogonal_symplectic(nmodes)
        @test isapprox(O', inv(O), atol = 1e-5)
        Omega = symplecticform(nmodes)
        @test isapprox(O' * Omega * O, Omega, atol = 1e-5)

        Spassive = randsymplectic(nmodes, passive = true)
        @test isapprox(Spassive', inv(Spassive), atol = 1e-5)
        @test isapprox(Spassive' * Omega * Spassive, Omega, atol = 1e-5)

        S = randsymplectic(nmodes)
        @test isapprox(S' * Omega * S, Omega, atol = 1e-5)
    end

    @testset "random states" begin
        nmodes = rand(1:20)
        rs = randstate(nmodes)
        rc = randchannel(nmodes)
        @test rc isa GaussianChannel
        @test rc * rs isa GaussianState
        Omega = symplecticform(nmodes)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rs.covar .+ im*Omega)))

        rspure = randstate(nmodes, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rspure.covar .+ im*Omega)))
        @test isapprox(purity(rspure), 1.0, atol = 1e-5)

        rs_float32 = randstate(AbstractArray{Float32}, nmodes)
        rc_float32 = randchannel(AbstractArray{Float32}, nmodes)
        @test rc_float32 isa GaussianChannel
        @test rc_float32 * rs_float32 isa GaussianState
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rs_float32.covar .+ im*Omega)))

        rspure_float32 = randstate(AbstractArray{Float32}, nmodes, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rspure_float32.covar .+ im*Omega)))
        @test isapprox(purity(rspure_float32), 1.0, atol = 1e-5)

        rs_static = randstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, nmodes)
        rc_static = randchannel(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, nmodes)
        @test rc_static isa GaussianChannel
        @test rc_static * rs_static isa GaussianState
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(Array(rs_static.covar) .+ im*Omega)))

        rspure_static = randstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, nmodes, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(Array(rspure_static.covar) .+ im*Omega)))
        @test isapprox(purity(rspure_float32), 1.0, atol = 1e-5)
    end

    @testset "random unitaries" begin
        nmodes = rand(1:20)
        ru = randunitary(nmodes)
        Omega = symplecticform(nmodes)
        @test isapprox(ru.symplectic' * Omega * ru.symplectic, Omega, atol = 1e-5)

        rupassive = randunitary(nmodes, passive = true)
        @test isapprox(rupassive.symplectic', inv(rupassive.symplectic), atol = 1e-5)
        @test isapprox(rupassive.symplectic' * Omega * rupassive.symplectic, Omega, atol = 1e-5)

        ru_float32 = randunitary(AbstractArray{Float32}, nmodes)
        @test isapprox(ru_float32.symplectic' * Omega * ru_float32.symplectic, Omega, atol = 1e-5)

        rupassive_float32 = randunitary(nmodes, passive = true)
        @test isapprox(rupassive_float32.symplectic', inv(rupassive_float32.symplectic), atol = 1e-5)
        @test isapprox(rupassive_float32.symplectic' * Omega * rupassive_float32.symplectic, Omega, atol = 1e-5)

        ru_static = randunitary(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, nmodes)
        @test isapprox(ru_static.symplectic' * Omega * ru_static.symplectic, Omega, atol = 1e-5)

        rupassive_static = randunitary(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, nmodes, passive = true)
        @test isapprox(rupassive_static.symplectic', inv(rupassive_static.symplectic), atol = 1e-5)
        @test isapprox(rupassive_static.symplectic' * Omega * rupassive_static.symplectic, Omega, atol = 1e-5)
    end

    @testset "random channels" begin
        #=
        nmodes = rand(1:20)
        rc = randchannel(nmodes)
        Omega = symplecticform(nmodes)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rc.noise .+ im*Omega .- im*rc.transform*Omega*rc.transform')))

        rc_float32 = randchannel(AbstractArray{Float32}, nmodes)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rc_float32.noise .+ im*Omega .- im*rc_float32.transform'*Omega*rc_float32.transform)))

        rc_static = randchannel(SVector{2*nmodes}, SMatrix{2*nmodes, 2*nmodes}, nmodes)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(Array(rc_static.noise) .+ im*Omega .- im*Array(rc_static.transform)'*Omega*Array(rc_static.transform))))
        =#
    end
end
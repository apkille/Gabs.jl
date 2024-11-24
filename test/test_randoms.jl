@testitem "CanonicalForm" begin
    using Gabs
    using StaticArrays
    using LinearAlgebra: eigvals, adjoint

    @testset "random utils" begin
        nmodes = rand(1:20)
        repr = CanonicalForm(nmodes)
        U = Gabs._rand_unitary(repr)
        @test isapprox(adjoint(U), inv(U), atol = 1e-5)

        O = Gabs._rand_orthogonal_symplectic(repr)
        @test isapprox(O', inv(O), atol = 1e-5)
        Omega = symplecticform(repr)
        @test isapprox(O' * Omega * O, Omega, atol = 1e-5)

        Spassive = randsymplectic(repr, passive = true)
        @test isapprox(Spassive', inv(Spassive), atol = 1e-5)
        @test isapprox(Spassive' * Omega * Spassive, Omega, atol = 1e-5)

        S = randsymplectic(repr)
        @test isapprox(S' * Omega * S, Omega, atol = 1e-5)

        S_array = randsymplectic(Array, repr)
        @test isapprox(S_array' * Omega * S_array, Omega, atol = 1e-5)
    end

    @testset "random states" begin
        nmodes = rand(1:20)
        repr = CanonicalForm(nmodes)
        rs = randstate(repr)
        rc = randchannel(repr)
        @test rc isa GaussianChannel
        @test rc * rs isa GaussianState
        Omega = symplecticform(repr)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rs.covar .+ im*Omega)))

        rspure = randstate(repr, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rspure.covar .+ im*Omega)))
        @test isapprox(purity(rspure), 1.0, atol = 1e-5)

        rs_array = randstate(Array, repr)
        rc_array = randchannel(Array, repr)
        @test rc_array isa GaussianChannel
        @test rc_array * rs_array isa GaussianState
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rs_array.covar .+ im*Omega)))

        rspure_array = randstate(Array, repr, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-3)), real(eigvals(rspure_array.covar .+ im*Omega)))
        @test isapprox(purity(rspure_array), 1.0, atol = 1e-3)

        rs_static = randstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, repr)
        rc_static = randchannel(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, repr)
        @test rc_static isa GaussianChannel
        @test rc_static * rs_static isa GaussianState
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(Array(rs_static.covar) .+ im*Omega)))

        rspure_static = randstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, repr, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(Array(rspure_static.covar) .+ im*Omega)))
        @test isapprox(purity(rspure_static), 1.0, atol = 1e-5)
    end

    @testset "random unitaries" begin
        nmodes = rand(1:20)
        repr = CanonicalForm(nmodes)
        ru = randunitary(repr)
        Omega = symplecticform(repr)
        @test isapprox(ru.symplectic' * Omega * ru.symplectic, Omega, atol = 1e-5)

        rupassive = randunitary(repr, passive = true)
        @test isapprox(rupassive.symplectic', inv(rupassive.symplectic), atol = 1e-5)
        @test isapprox(rupassive.symplectic' * Omega * rupassive.symplectic, Omega, atol = 1e-5)

        ru_array = randunitary(Array, repr)
        @test isapprox(ru_array.symplectic' * Omega * ru_array.symplectic, Omega, atol = 1e-5)

        rupassive_array = randunitary(repr, passive = true)
        @test isapprox(rupassive_array.symplectic', inv(rupassive_array.symplectic), atol = 1e-5)
        @test isapprox(rupassive_array.symplectic' * Omega * rupassive_array.symplectic, Omega, atol = 1e-5)

        ru_static = randunitary(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, repr)
        @test isapprox(ru_static.symplectic' * Omega * ru_static.symplectic, Omega, atol = 1e-5)

        rupassive_static = randunitary(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, repr, passive = true)
        @test isapprox(rupassive_static.symplectic', inv(rupassive_static.symplectic), atol = 1e-5)
        @test isapprox(rupassive_static.symplectic' * Omega * rupassive_static.symplectic, Omega, atol = 1e-5)
    end

    @testset "random channels" begin
        nmodes = rand(1:20)
        repr = CanonicalForm(nmodes)
        rc = randchannel(repr)
        Omega = symplecticform(repr)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rc.noise .+ im*Omega .- im*rc.transform*Omega*rc.transform')))

        rc_array = randchannel(Array, repr)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rc_array.noise .+ im*Omega .- im*rc_array.transform*Omega*rc_array.transform')))

        rc_static = randchannel(SVector{2*nmodes}, SMatrix{2*nmodes, 2*nmodes}, repr)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(Array(rc_static.noise) .+ im*Omega .- im*Array(rc_static.transform)*Omega*Array(rc_static.transform)')))
    end
end
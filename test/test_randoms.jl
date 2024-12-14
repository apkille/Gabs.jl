@testitem "Random objects" begin
    using Gabs
    using StaticArrays
    using LinearAlgebra: eigvals, adjoint

    @testset "random utils" begin
        nmodes = rand(1:20)
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)
        U_qpair = Gabs._rand_unitary(qpairbasis)
        U_qblock = Gabs._rand_unitary(qblockbasis)
        @test isapprox(adjoint(U_qpair), inv(U_qpair), atol = 1e-5)
        @test isapprox(adjoint(U_qblock), inv(U_qblock), atol = 1e-5)

        O_qpair = Gabs._rand_orthogonal_symplectic(qpairbasis)
        O_qblock = Gabs._rand_orthogonal_symplectic(qblockbasis)
        @test isapprox(O_qpair', inv(O_qpair), atol = 1e-5)
        @test isapprox(O_qblock', inv(O_qblock), atol = 1e-5)
        Omega_qpair = symplecticform(qpairbasis)
        Omega_qblock= symplecticform(qblockbasis)
        @test isapprox(O_qpair' * Omega_qpair * O_qpair, Omega_qpair, atol = 1e-5)
        @test isapprox(O_qblock' * Omega_qblock * O_qblock, Omega_qblock, atol = 1e-5)

        Spassive_qpair = randsymplectic(qpairbasis, passive = true)
        Spassive_qblock = randsymplectic(qblockbasis, passive = true)
        @test isapprox(Spassive_qpair', inv(Spassive_qpair), atol = 1e-5)
        @test isapprox(Spassive_qpair' * Omega_qpair * Spassive_qpair, Omega_qpair, atol = 1e-5)
        @test isapprox(Spassive_qblock', inv(Spassive_qblock), atol = 1e-5)
        @test isapprox(Spassive_qblock' * Omega_qblock * Spassive_qblock, Omega_qblock, atol = 1e-5)

        S_qpair = randsymplectic(qpairbasis)
        S_qblock = randsymplectic(qblockbasis)
        @test isapprox(S_qpair' * Omega_qpair * S_qpair, Omega_qpair, atol = 1e-5)
        @test isapprox(S_qblock' * Omega_qblock * S_qblock, Omega_qblock, atol = 1e-5)

        S_array = randsymplectic(Array, qpairbasis)
        @test isapprox(S_array' * Omega_qpair * S_array, Omega_qpair, atol = 1e-5)
    end

    @testset "random states" begin
        nmodes = rand(1:20)
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)
        rs = randstate(qpairbasis)
        rc = randchannel(qpairbasis)
        @test rc isa GaussianChannel
        @test rc * rs isa GaussianState
        Omega = symplecticform(qpairbasis)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rs.covar .+ im*Omega)))

        rspure = randstate(qpairbasis, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rspure.covar .+ im*Omega)))
        @test isapprox(purity(rspure), 1.0, atol = 1e-5)

        rs_array = randstate(Array, qpairbasis)
        rc_array = randchannel(Array, qpairbasis)
        @test rc_array isa GaussianChannel
        @test rc_array * rs_array isa GaussianState
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rs_array.covar .+ im*Omega)))

        rspure_array = randstate(Array, qpairbasis, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-3)), real(eigvals(rspure_array.covar .+ im*Omega)))
        @test isapprox(purity(rspure_array), 1.0, atol = 1e-3)

        rs_static = randstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis)
        rc_static = randchannel(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis)
        @test rc_static isa GaussianChannel
        @test rc_static * rs_static isa GaussianState
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(Array(rs_static.covar) .+ im*Omega)))

        rspure_static = randstate(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, pure = true)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(Array(rspure_static.covar) .+ im*Omega)))
        @test isapprox(purity(rspure_static), 1.0, atol = 1e-5)
    end

    @testset "random unitaries" begin
        nmodes = rand(1:20)
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)
        ru = randunitary(qpairbasis)
        Omega = symplecticform(qpairbasis)
        @test isapprox(ru.symplectic' * Omega * ru.symplectic, Omega, atol = 1e-5)

        rupassive = randunitary(qpairbasis, passive = true)
        @test isapprox(rupassive.symplectic', inv(rupassive.symplectic), atol = 1e-5)
        @test isapprox(rupassive.symplectic' * Omega * rupassive.symplectic, Omega, atol = 1e-5)

        ru_array = randunitary(Array, qpairbasis)
        @test isapprox(ru_array.symplectic' * Omega * ru_array.symplectic, Omega, atol = 1e-5)

        rupassive_array = randunitary(qpairbasis, passive = true)
        @test isapprox(rupassive_array.symplectic', inv(rupassive_array.symplectic), atol = 1e-5)
        @test isapprox(rupassive_array.symplectic' * Omega * rupassive_array.symplectic, Omega, atol = 1e-5)

        ru_static = randunitary(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis)
        @test isapprox(ru_static.symplectic' * Omega * ru_static.symplectic, Omega, atol = 1e-5)

        rupassive_static = randunitary(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, passive = true)
        @test isapprox(rupassive_static.symplectic', inv(rupassive_static.symplectic), atol = 1e-5)
        @test isapprox(rupassive_static.symplectic' * Omega * rupassive_static.symplectic, Omega, atol = 1e-5)
    end

    @testset "random channels" begin
        nmodes = rand(1:20)
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)
        rc = randchannel(qpairbasis)
        Omega = symplecticform(qpairbasis)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rc.noise .+ im*Omega .- im*rc.transform*Omega*rc.transform')))

        rc_array = randchannel(Array, qpairbasis)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(rc_array.noise .+ im*Omega .- im*rc_array.transform*Omega*rc_array.transform')))

        rc_static = randchannel(SVector{2*nmodes}, SMatrix{2*nmodes, 2*nmodes}, qpairbasis)
        @test all(i -> ((i >= 0) || isapprox(i, 0.0, atol = 1e-5)), real(eigvals(Array(rc_static.noise) .+ im*Omega .- im*Array(rc_static.transform)*Omega*Array(rc_static.transform)')))
    end
end
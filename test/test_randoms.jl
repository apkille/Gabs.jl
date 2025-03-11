@testitem "Random objects" begin
    using Gabs
    using StaticArrays
    using LinearAlgebra: eigvals, adjoint

    @testset "random utils" begin
        nmodes = rand(1:5)
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
        @test issymplectic(qpairbasis, O_qpair, atol = 1e-5)
        @test issymplectic(qblockbasis, O_qblock, atol = 1e-5)

        Spassive_qpair = randsymplectic(qpairbasis, passive = true)
        Spassive_qblock = randsymplectic(qblockbasis, passive = true)
        @test isapprox(Spassive_qpair', inv(Spassive_qpair), atol = 1e-5)
        @test issymplectic(qpairbasis, Spassive_qpair, atol = 1e-5)
        @test isapprox(Spassive_qblock', inv(Spassive_qblock), atol = 1e-5)
        @test issymplectic(qblockbasis, Spassive_qblock, atol = 1e-5)

        S_qpair = randsymplectic(qpairbasis)
        S_qblock = randsymplectic(qblockbasis)
        @test issymplectic(qpairbasis, S_qpair, atol = 1e-5)
        @test issymplectic(qblockbasis, S_qblock, atol = 1e-5)

        S_array = randsymplectic(Array, qpairbasis)
        @test issymplectic(qpairbasis, S_array, atol = 1e-5)
    end

    @testset "random states" begin
        nmodes, ħ = rand(1:5), rand(1:5)
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)
        rs_pair = randstate(qpairbasis)
        rs_block = randstate(qblockbasis)
        rc = randchannel(qpairbasis)
        @test rc isa GaussianChannel
        @test rc * rs_pair isa GaussianState
        @test rs_pair.ħ == 2 && rs_block.ħ == 2
        @test isgaussian(rs_pair, atol = 1e-5)
        @test isgaussian(rs_block, atol = 1e-5)

        rspure_pair = randstate(qpairbasis, pure = true, ħ = ħ)
        rspure_block = randstate(qblockbasis, pure = true, ħ = ħ)
        @test isgaussian(rspure_pair, atol = 1e-5)
        @test isapprox(purity(rspure_pair), 1.0, atol = 1e-5)
        @test isgaussian(rspure_block, atol = 1e-5)
        @test isapprox(purity(rspure_block), 1.0, atol = 1e-5)

        rs_array, rs_array_static  = randstate(Array, qpairbasis), randstate(SArray, qpairbasis)
        rc_array, rc_array_static = randchannel(Array, qpairbasis),  randchannel(SArray, qpairbasis)
        @test rc_array isa GaussianChannel &&  rc_array_static isa GaussianChannel
        @test rc_array.ħ == 2 && rc_array_static.ħ == 2
        @test rc_array * rs_array isa GaussianState && rc_array_static * rs_array_static isa GaussianState
        @test isgaussian(rs_array, atol = 1e-5) && isgaussian(rs_array_static, atol = 1e-5)

        rspure_array = randstate(SArray, qpairbasis, pure = true)
        @test isgaussian(rspure_array, atol = 1e-5)
        @test_broken isapprox(purity(rspure_array), 1.0, atol = 1e-3)

        rs_static, rs_static1, rs_static2,  = randstate(SVector{2*nmodes, Float64}, SMatrix{2*nmodes,2*nmodes, Float64}, qpairbasis), randstate(SVector, SMatrix, qpairbasis), randstate(SArray, qpairbasis)
        rc_static, rc_static1, rc_static2 = randchannel(SVector{2*nmodes, Float64}, SMatrix{2*nmodes,2*nmodes, Float64}, qpairbasis), randchannel(SVector, SMatrix, qpairbasis), randchannel(SArray, qpairbasis)
        @test rc_static isa GaussianChannel && rc_static1 isa GaussianChannel && rc_static2 isa GaussianChannel
        @test rs_static isa GaussianState && rs_static1 isa GaussianState && rs_static2 isa GaussianState
        @test rc_static * rs_static isa GaussianState && rc_static1 * rs_static1 isa GaussianState && rc_static2 * rs_static2 isa GaussianState
        @test isgaussian(rs_static, atol = 1e-5) && isgaussian(rs_static1, atol = 1e-5) && isgaussian(rs_static2, atol = 1e-5)

        rspure_static = randstate(SVector{2*nmodes, Float64}, SMatrix{2*nmodes,2*nmodes, Float64}, qpairbasis, pure = true)
        @test isgaussian(rspure_static, atol = 1e-5)
        @test isapprox(purity(rspure_static), 1.0, atol = 1e-5)
    end

    @testset "random unitaries" begin
        nmodes, ħ = rand(1:5), rand(1:5)
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)
        ru = randunitary(qpairbasis, ħ = ħ)
        @test isgaussian(ru, atol = 1e-5)

        rupassive = randunitary(qpairbasis, passive = true)
        @test rupassive.ħ == 2
        @test isapprox(rupassive.symplectic', inv(rupassive.symplectic), atol = 1e-5)
        @test isgaussian(rupassive, atol = 1e-5)

        ru_array, ru_array_static = randunitary(Array, qpairbasis), randunitary(SArray, qpairbasis)
        @test ru_array.ħ == 2 && ru_array_static.ħ == 2
        @test isgaussian(ru_array, atol = 1e-5) && isgaussian(ru_array_static, atol = 1e-5)

        rupassive_array = randunitary(qpairbasis, passive = true)
        @test isapprox(rupassive_array.symplectic', inv(rupassive_array.symplectic), atol = 1e-5)
        @test isgaussian(rupassive_array, atol = 1e-5)

        ru_static = randunitary(SVector, SMatrix, qpairbasis)
        @test ru_static.ħ == 2
        @test isgaussian(ru_static, atol = 1e-5)

        rupassive_static = randunitary(SVector{2*nmodes, Float64}, SMatrix{2*nmodes, 2*nmodes, Float64}, qpairbasis, passive = true)
        @test isapprox(rupassive_static.symplectic', inv(rupassive_static.symplectic), atol = 1e-5)
        @test isgaussian(rupassive_static, atol = 1e-5)
    end

    @testset "random channels" begin
        nmodes, ħ = rand(1:5), rand(1:5)
        qpairbasis = QuadPairBasis(nmodes)
        qblockbasis = QuadBlockBasis(nmodes)
        rc = randchannel(qpairbasis, ħ = ħ)
        @test isgaussian(rc, atol = 1e-5)

        rc_array, rc_array_static = randchannel(Array, qpairbasis), randchannel(SArray, qpairbasis)
        @test rc_array.ħ == 2 && rc_array_static.ħ == 2
        @test isgaussian(rc_array, atol = 1e-5) && isgaussian(rc_array_static, atol = 1e-5)

        rc_static = randchannel(SVector, SMatrix, qpairbasis)
        @test rc_static.ħ == 2
        @test isgaussian(rc_static, atol = 1e-5)
    end
end
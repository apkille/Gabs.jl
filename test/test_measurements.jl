@testitem "Measurements" begin
    using Gabs
    using StaticArrays
    using LinearAlgebra: det, I

    @testset "generaldyne" begin
        qpairbasis, qblockbasis = QuadPairBasis(1), QuadBlockBasis(1)
        for basis in [qpairbasis, qblockbasis]
            vac = vacuumstate(basis)
            vacs = vac ⊗ vac ⊗ vac ⊗ vac
            gd1 = generaldyne(vacs, [2, 4], proj = vac ⊗ vac)
            @test isapprox(gd1.result, vac ⊗ vac, atol = 1e-12)
            @test isapprox(gd1.state, vacs, atol = 1e-12)

            coh = coherentstate(basis, 1.0+im)
            cohs = coh ⊗ vac ⊗ coh ⊗ vac
            epr = eprstate(basis ⊕ basis, 1.0, 3.0)
            gd2 = generaldyne(cohs, [1, 4], proj = epr)
            @test isapprox(gd2.result, epr, atol = 1e-12)
            @test isapprox(gd2.state, vac ⊗ vac ⊗ coh ⊗ vac, atol = 1e-12)
        end

        indices, nmodes = [7, 8, 9, 10], 10

        # random Gaussian state tests
        rs_qpair = randstate(QuadPairBasis(nmodes))
        rs_qblock = changebasis(QuadBlockBasis, rs_qpair)
        proj_qpair = randstate(QuadPairBasis(length(indices)))
        proj_qblock = changebasis(QuadBlockBasis, proj_qpair)
        gd3_qpair = generaldyne(rs_qpair, indices, proj = proj_qpair)
        gd3_qblock = generaldyne(rs_qblock, indices, proj = proj_qblock)

        # analytical calculation of evolved subsystem that is not measured
        xA, xB, VA, VB, VAB = Gabs._part_state(rs_qpair, indices)
        gd3_evolved_mean = xA .+ VAB * ((inv(VB .+ proj_qpair.covar)) * (proj_qpair.mean .- xB))
        gd3_evolved_covar = VA .- VAB * ((inv(VB .+ proj_qpair.covar)) * transpose(VAB))
        out3_mean = vcat(gd3_evolved_mean, zeros(2*length(indices)))
        out3_covar = zeros(2*nmodes, 2*nmodes)
        copyto!(@view(out3_covar[1:2*(nmodes-length(indices)), 1:2*(nmodes-length(indices))]), gd3_evolved_covar)
        copyto!(@view(out3_covar[2*(nmodes-length(indices))+1:2*nmodes, 2*(nmodes-length(indices))+1:2*nmodes]), Matrix{Float64}(I, 2*length(indices), 2*length(indices)))

        evolved_state_qpair = GaussianState(QuadPairBasis(nmodes), out3_mean, out3_covar)
        evolved_state_qblock = changebasis(QuadBlockBasis, evolved_state_qpair)
        @test isapprox(gd3_qpair.state, evolved_state_qpair)
        @test isapprox(gd3_qblock.state, evolved_state_qblock)
        
        # tests to check that static arrays are outputted for generaldyne 
        # measurements of Gaussian states wrapping static arrays
        sstatic = vacuumstate(SVector{2}, SMatrix{2,2}, QuadPairBasis(1))
        statestatic = sstatic ⊗ sstatic ⊗ sstatic ⊗ sstatic
        gdstatic = generaldyne(statestatic, [2])
        @test (gdstatic.state).mean isa SVector && (gdstatic.state).covar isa SMatrix
        @test isequal(gdstatic.state, statestatic)
    end
end
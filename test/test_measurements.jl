@testitem "Measurements" begin
    using Gabs
    using StaticArrays
    using LinearAlgebra: det

    @testset "generaldyne" begin
        basis = QuadPairBasis(1)
        vac = vacuumstate(basis)
        vacs = vac ⊗ vac ⊗ vac ⊗ vac
        gd1 = Generaldyne(vacs, vac ⊗ vac, [2, 4])
        out1 = output(gd1)
        @test isequal(out1, vac ⊗ vac)

        coh = coherentstate(basis, 1.0+im)
        cohs = coh ⊗ vac ⊗ coh ⊗ vac
        epr = eprstate(basis ⊕ basis, 1.0, 3.0)
        gd2 = Generaldyne(cohs, epr, [1, 4])
        out2 = output(gd2)
        @test isequal(out2, vac ⊗ coh)
        out2_prob = exp(transpose(epr.mean .- (coh ⊗ vac).mean) * (inv(epr.covar .+ (coh ⊗ vac).covar) * (epr.mean .- (coh ⊗ vac).mean)))/(pi^2 * sqrt(det(epr.covar .+ (coh ⊗ vac).covar)))
        @test isapprox(prob(gd2), out2_prob)

        state = GaussianState(2*basis, Vector{Float64}(collect(1:4)), Matrix{Float64}(reshape(collect(1:16), (4,4))))
        meas = GaussianState(basis, Vector{Float64}(collect(1:2)), Matrix{Float64}(reshape(collect(1:4), (2,2))))
        gd3 = Generaldyne(state, meas, [2])
        out3 = output(gd3)
        xA, xB = [1.0, 2.0], [3.0, 4.0]
        VA, VB, VAB = [1.0 5.0; 2.0 6.0], [11.0 15.0; 12.0 16.0], [9.0 13.0; 10.0 14.0]
        out3_mean = xA .+ VAB*((inv(VB .+ meas.covar))*(meas.mean .- xB))
        out3_covar = VA .- VAB*((inv(VB .+ meas.covar))*transpose(VAB))
        @test isapprox(out3, GaussianState(basis, out3_mean, out3_covar))

        sstatic = vacuumstate(SVector{2}, SMatrix{2,2}, basis)
        statestatic = sstatic ⊗ sstatic ⊗ sstatic ⊗ sstatic
        gdstatic = Generaldyne(statestatic, sstatic, [2])
        outstatic = output(gdstatic)
        @test isequal(outstatic, sstatic ⊗ sstatic ⊗ sstatic)
    end
end
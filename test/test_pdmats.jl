@testitem "PDMats Integration" begin
    using Gabs
    using BlockArrays
    using LinearAlgebra: norm, det
    using PDMats
    
    @testset "PDGaussianState creation" begin
        basis = QuadPairBasis(2)
        
        std_state = vacuumstate(basis)
        pd_state = vacuumstate_pd(basis)
        
        @test isapprox(Matrix(pd_state.covar), std_state.covar)
        @test isapprox(pd_state.mean, std_state.mean)
        
        thermal_std = thermalstate(basis, 2.0)
        thermal_pd = thermalstate_pd(basis, 2.0)
        
        @test isapprox(Matrix(thermal_pd.covar), thermal_std.covar)
        @test isapprox(thermal_pd.mean, thermal_std.mean)
        
        alpha = 1.0 + im
        coherent_std = coherentstate(basis, alpha)
        coherent_pd = coherentstate_pd(basis, alpha)
        
        @test isapprox(Matrix(coherent_pd.covar), coherent_std.covar)
        @test isapprox(coherent_pd.mean, coherent_std.mean)
    end
    
    @testset "Core operations" begin
        basis = QuadPairBasis(3)
        std_state = thermalstate(basis, 1.5)
        pd_state = thermalstate_pd(basis, 1.5)
        testvec = randn(2*basis.nmodes)
        
        @test isapprox(purity(std_state), purity(pd_state))
        
        std_whitened = whiten(std_state, testvec)
        pd_whitened = whiten(pd_state, testvec)
        @test isapprox(std_whitened, pd_whitened)
        
        std_unwhitened = unwhiten(std_state, testvec)
        pd_unwhitened = unwhiten(pd_state, testvec)
        @test isapprox(std_unwhitened, pd_unwhitened)
        
        @test isapprox(quad(std_state, testvec), quad(pd_state, testvec))
        @test isapprox(invquad(std_state, testvec), invquad(pd_state, testvec))
    end
    
    @testset "Constructor from GaussianState" begin
        basis = QuadPairBasis(2)
        
        rand_state = randstate(basis)
        
        pd_state = PDGaussianState(rand_state)
        
        @test isapprox(Matrix(pd_state.covar), rand_state.covar)
        @test isapprox(pd_state.mean, rand_state.mean)
        @test pd_state.basis == rand_state.basis
        @test pd_state.ħ == rand_state.ħ
    end
end
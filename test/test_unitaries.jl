@testitem "Unitaries" begin
    using Gabs
    using StaticArrays
    
    @testset "displacement operator" begin
        alpha = rand(ComplexF64)
        @test displace(alpha) isa GaussianUnitary
        @test displace(SVector{2}, SMatrix{2,2}, alpha) isa GaussianUnitary
    end

    @testset "squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        @test squeeze(r, theta) isa GaussianUnitary
        @test squeeze(SVector{2}, SMatrix{2,2}, r, theta) isa GaussianUnitary
    end

    @testset "two-mode squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        @test twosqueeze(r, theta) isa GaussianUnitary
        @test twosqueeze(SVector{4}, SMatrix{4,4}, r, theta) isa GaussianUnitary
    end

    @testset "phase-shift operator" begin
        theta = rand(Float64)
        @test phaseshift(theta) isa GaussianUnitary
        @test phaseshift(SVector{2}, SMatrix{2,2}, theta) isa GaussianUnitary
    end

    @testset "beamsplitter operator" begin
        theta = rand(Float64)
        @test beamsplitter(theta) isa GaussianUnitary
        @test beamsplitter(SVector{4}, SMatrix{4,4}, theta) isa GaussianUnitary
    end

    @testset "tensor products" begin
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(alpha1), displace(alpha2)
        ds = tensor(d1, d2)
        @test ds isa GaussianUnitary
        @test ds == d1 ⊗ d2
        @test tensor(SVector{4}, SMatrix{4,4}, d1, d2) isa GaussianUnitary

        r, theta = rand(Float64), rand(Float64)
        p = phaseshift(theta)
        @test tensor(p, tensor(d1, d2)) == p ⊗ d1 ⊗ d2


        dstatic = displace(SVector{2}, SMatrix{2,2}, alpha1)
        tpstatic = dstatic ⊗ dstatic ⊗ dstatic
        @test tpstatic.disp isa SVector{6}
        @test tpstatic.symplectic isa SMatrix{6,6}
        tp = dstatic ⊗ d1 ⊗ dstatic
        @test tp.disp isa Vector
        @test tp.symplectic isa Matrix
    end

    @testset "actions" begin
        alpha = rand(ComplexF64)
        d = displace(alpha)
        v = vacuumstate()
        c = coherentstate(alpha)
        @test d * v == c
        @test apply!(v, d) == c
        
        v1, v2 = vacuumstate(), vacuumstate()
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(alpha1), displace(alpha2)
        c1, c2 = coherentstate(alpha1), coherentstate(alpha2)
        @test (d1 ⊗ d2) * (v1 ⊗ v2) == c1 ⊗ c2
    end
end
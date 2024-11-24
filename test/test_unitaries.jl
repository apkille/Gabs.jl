@testitem "CanonicalForm" begin
    using Gabs
    using StaticArrays
    
    repr1 = CanonicalForm(1)

    @testset "displacement operator" begin
        alpha = rand(ComplexF64)
        @test displace(repr1, alpha) isa GaussianUnitary
        @test displace(Array, repr1, alpha) isa GaussianUnitary
        @test displace(SVector{2}, SMatrix{2,2}, repr1, alpha) isa GaussianUnitary
    end

    @testset "squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        @test squeeze(repr1, r, theta) isa GaussianUnitary
        @test squeeze(Array, repr1, r, theta) isa GaussianUnitary
        @test squeeze(SVector{2}, SMatrix{2,2}, repr1, r, theta) isa GaussianUnitary
    end

    @testset "two-mode squeeze operator" begin
        r, theta = rand(Float64), rand(Float64)
        @test twosqueeze(2*repr1, r, theta) isa GaussianUnitary
        @test twosqueeze(Array, 2*repr1, r, theta) isa GaussianUnitary
        @test twosqueeze(SVector{4}, SMatrix{4,4}, 2*repr1, r, theta) isa GaussianUnitary
    end

    @testset "phase-shift operator" begin
        theta = rand(Float64)
        @test phaseshift(repr1, theta) isa GaussianUnitary
        @test phaseshift(Array, repr1, theta) isa GaussianUnitary
        @test phaseshift(SVector{2}, SMatrix{2,2}, repr1, theta) isa GaussianUnitary
    end

    @testset "beamsplitter operator" begin
        theta = rand(Float64)
        @test beamsplitter(2*repr1, theta) isa GaussianUnitary
        @test beamsplitter(Array, 2*repr1, theta) isa GaussianUnitary
        @test beamsplitter(SVector{4}, SMatrix{4,4}, 2*repr1, theta) isa GaussianUnitary
    end

    @testset "tensor products" begin
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(repr1, alpha1), displace(repr1, alpha2)
        ds = tensor(d1, d2)
        @test ds isa GaussianUnitary
        @test ds == d1 ⊗ d2
        @test isapprox(ds, d1 ⊗ d2)
        @test tensor(SVector{4}, SMatrix{4,4}, d1, d2) isa GaussianUnitary

        r, theta = rand(Float64), rand(Float64)
        p = phaseshift(repr1, theta)
        @test tensor(p, tensor(d1, d2)) == p ⊗ d1 ⊗ d2

        dstatic = displace(SVector{2}, SMatrix{2,2}, repr1, alpha1)
        tpstatic = dstatic ⊗ dstatic ⊗ dstatic
        @test tpstatic.disp isa SVector{6}
        @test tpstatic.symplectic isa SMatrix{6,6}
        tp = dstatic ⊗ d1 ⊗ dstatic
        @test tp.disp isa Vector
        @test tp.symplectic isa Matrix
    end

    @testset "actions" begin
        alpha = rand(ComplexF64)
        d = displace(repr1, alpha)
        v = vacuumstate(repr1)
        c = coherentstate(repr1, alpha)
        @test d * v == c
        @test apply!(v, d) == c
        
        v1, v2 = vacuumstate(repr1), vacuumstate(repr1)
        alpha1, alpha2 = rand(ComplexF64), rand(ComplexF64)
        d1, d2 = displace(repr1, alpha1), displace(repr1, alpha2)
        c1, c2 = coherentstate(repr1, alpha1), coherentstate(repr1, alpha2)
        @test (d1 ⊗ d2) * (v1 ⊗ v2) == c1 ⊗ c2
    end
end
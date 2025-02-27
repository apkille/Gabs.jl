@testitem "Symbolic Channels" begin
    using Gabs
    using Symbolics
    using StaticArrays

    nmodes = 2*rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)
    noise = zeros(Num, 2*nmodes, 2*nmodes)

    @variables a b
    @variables as[1:nmodes] bs[1:nmodes]
    as, bs = collect(as), collect(bs)

    test_configs = [
    (displace, (a, noise), (as, noise), "symbolic displacement operator"), 
    (squeeze, (a, b, noise), (as, bs, noise), "symbolic squeeze operator"), 
    (phaseshift, (a, noise), (as, noise), "symbolic phase-shift operator"),
    (twosqueeze, (a, b, noise), (as, bs, noise), "symbolic two-mode squeeze operator"),
    (beamsplitter, (a, noise), (as, noise), "symbolic beamsplitter operator"),
    (amplifier, (a, b), (as, bs), "symbolic amplifier channel"),
    (attenuator, (a, b), (as, bs), "symbolic attenuator channel")
    ]

    for (op, params, multi_params, desc) in test_configs
        @testset "$desc" begin
            op_pair = op(qpairbasis, params...)
            op_block = op(qblockbasis, params...)  
            @test op_pair isa GaussianChannel && op_block isa GaussianChannel
            @test op(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, params...) isa GaussianChannel
            @test op(Array, qpairbasis, params...) isa GaussianChannel
            @test isapprox(op_block, changebasis(QuadBlockBasis, op_pair))

            op_pair_multi = op(qpairbasis, multi_params...)
            op_block_multi = op(qblockbasis, multi_params...)
            @test op_pair_multi isa GaussianChannel && op_block_multi isa GaussianChannel
            @test op(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, multi_params...) isa GaussianChannel
            @test op(Array, qpairbasis, multi_params...) isa GaussianChannel
            @test isapprox(op_block_multi, changebasis(QuadBlockBasis, op_pair_multi))
        end
    end

    @testset "Symbolic tensor products" begin
        @variables alpha1 alpha2 theta
        noise_ds = [noise zeros(2*nmodes,2*nmodes); zeros(2*nmodes,2*nmodes) noise]
        d1 = displace(qpairbasis, alpha1, noise)
        d2 = displace(qpairbasis, alpha2, noise)
        ds = tensor(d1, d2)
        @test ds isa GaussianChannel
        @test isapprox(ds, d1 ⊗ d2)
        @test isapprox(ds, d1 ⊗ d2, atol = 1e-10)
        @test tensor(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, d1, d2) isa GaussianChannel
        p = phaseshift(qpairbasis, theta, noise)
        @test isapprox(tensor(tensor(p, d1), d2), p ⊗ d1 ⊗ d2)
        p_block = phaseshift(qblockbasis, theta, noise)
        p_blocks = phaseshift(2*qblockbasis, repeat([theta], 2*nmodes), noise_ds)
        @test isapprox(p_block ⊗ p_block, p_blocks)
        @variables alphas1[1:nmodes] alphas2[1:nmodes] thetas[1:nmodes]
        alphas1_vec = collect(alphas1)
        alphas2_vec = collect(alphas2)
        thetas_vec = collect(thetas)
        d1 = displace(qpairbasis, alphas1_vec, noise)
        d2 = displace(qpairbasis, alphas2_vec, noise)
        ds = tensor(d1, d2)
        @test ds isa GaussianChannel
        @test isapprox(ds, d1 ⊗ d2)
        @test isapprox(ds, d1 ⊗ d2, atol = 1e-10)
        @test tensor(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, d1, d2) isa GaussianChannel
        p = phaseshift(qpairbasis, thetas_vec, noise)
        @test isapprox(tensor(tensor(p, d1), d2), p ⊗ d1 ⊗ d2)
    end

    @testset "Symbolic actions" begin
        @variables alpha1 alpha2
        d = displace(qpairbasis, alpha1, noise)
        v = vacuumstate(qpairbasis)
        c = coherentstate(qpairbasis, alpha1)
        @test isapprox(d * v, c)
        @test isapprox(d * v, c, atol = 1e-10)
        @test_broken isapprox(apply!(v, d), c)
        v1, v2 = vacuumstate(qpairbasis), vacuumstate(qpairbasis)
        d1 = displace(qpairbasis, alpha1, noise)
        d2 = displace(qpairbasis, alpha2, noise)
        c1 = coherentstate(qpairbasis, alpha1)
        c2 = coherentstate(qpairbasis, alpha2)
        @test isapprox((d1 ⊗ d2) * (v1 ⊗ v2), c1 ⊗ c2)
    end
end

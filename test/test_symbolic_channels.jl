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
end

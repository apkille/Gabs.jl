@testitem "Symbolic Channels" begin
    using Gabs
    using Symbolics
    using StaticArrays

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)
    noise = zeros(Num, 2*nmodes, 2*nmodes)
    noise_ds = [noise zeros(2*nmodes,2*nmodes); zeros(2*nmodes,2*nmodes) noise]

    @testset "displacement operator" begin
        @variables α
        @variables αs[1:2*nmodes]
        alphas_vec = collect(αs)
        op_pair = displace(qpairbasis, α, noise)
        op_block = displace(qblockbasis, α, noise)
        @test op_pair isa GaussianChannel && op_block isa GaussianChannel
        @test displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, α, noise) isa GaussianChannel
        @test displace(Array, qpairbasis, α, noise) isa GaussianChannel
        @test isapprox(displace(qblockbasis, α, noise), changebasis(QuadBlockBasis, op_pair))
        @test isapprox(displace(qblockbasis, alphas_vec, noise), changebasis(QuadBlockBasis, displace(qpairbasis, alphas_vec, noise)))
    end

    @testset "Symbolic Squeeze and Phase-Shift Operators" begin
        @variables r θs
        @variables rs[1:nmodes] thetas_s[1:nmodes]
        rs_vec = collect(rs)
        thetas_vec = collect(thetas_s)
        @variables θps
        @variables thetas_ps[1:nmodes]
        thetas_ps_vec = collect(thetas_ps)
        test_configs = [(squeeze, (r, θs), (rs_vec, thetas_vec), "Squeeze Operator"), (phaseshift, (θps,), (thetas_ps_vec,), "Phase-Shift Operator")]
        for (op, single_params, multi_params, desc) in test_configs
            @testset "$desc" begin
                op_pair   = op(qpairbasis, single_params..., noise)
                op_block  = op(qblockbasis, single_params..., noise)
                @test op_pair isa GaussianChannel && op_block isa GaussianChannel
                @test iszero(simplify(op_block.transform - changebasis(QuadBlockBasis, op_pair).transform))
                @test iszero(simplify(op_block.disp - changebasis(QuadBlockBasis, op_pair).disp))
                @test op(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, single_params..., noise) isa GaussianChannel
                @test op(Array, qpairbasis, single_params..., noise) isa GaussianChannel
                @test isapprox(op(qblockbasis, single_params..., noise), changebasis(QuadBlockBasis, op_pair))
            end
            @testset "$desc" begin
                op_pair_multi  = op(qpairbasis, multi_params..., noise)
                op_block_multi = op(qblockbasis, multi_params..., noise)
                @test op_pair_multi isa GaussianChannel && op_block_multi isa GaussianChannel
                @test iszero(simplify(op_block_multi.transform - changebasis(QuadBlockBasis, op_pair_multi).transform))
                @test iszero(simplify(op_block_multi.disp - changebasis(QuadBlockBasis, op_pair_multi).disp))
                @test op(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, multi_params..., noise) isa GaussianChannel
                @test op(Array, qpairbasis, multi_params..., noise) isa GaussianChannel
                @test isapprox(op(qblockbasis, multi_params..., noise), changebasis(QuadBlockBasis, op(qpairbasis, multi_params..., noise)))
            end
        end
    end

    @testset "Symbolic Two-Mode Squeeze and Beamsplitter Operators" begin
        @variables r θ θ_bs
        @variables rs[1:nmodes] thetas[1:nmodes] thetas_bs[1:nmodes]
        rs_vec = collect(rs)
        thetas_vec = collect(thetas)
        thetas_bs_vec = collect(thetas_bs)
        test_configs = [(twosqueeze, (r, θ), (rs_vec, thetas_vec), "Two-Mode Squeeze Operator"),(beamsplitter, (θ_bs,), (thetas_bs_vec,), "Beamsplitter Operator")]
        for (op, single_params, multi_params, desc) in test_configs
            @testset "$desc" begin
                op_pair   = op(2*qpairbasis, single_params..., noise_ds)
                op_block  = op(2*qblockbasis, single_params..., noise_ds)
                @test op_pair isa GaussianChannel && op_block isa GaussianChannel
                @test isapprox(op_block, changebasis(QuadBlockBasis, op_pair))
                @test op(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, single_params..., noise_ds) isa GaussianChannel
                @test op(Array, 2*qpairbasis, single_params..., noise_ds) isa GaussianChannel
            end
            @testset "$desc" begin
                op_pair_multi  = op(2*qpairbasis, multi_params..., noise_ds)
                op_block_multi = op(2*qblockbasis, multi_params..., noise_ds)
                @test op_pair_multi isa GaussianChannel && op_block_multi isa GaussianChannel
                @test op(Array, 2*qpairbasis, multi_params..., noise_ds) isa GaussianChannel
                @test op(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, multi_params..., noise_ds) isa GaussianChannel
                @test isapprox(op_block_multi, changebasis(QuadBlockBasis, op_pair_multi))
            end
        end
    end

    @testset "Symbolic Gaussian Channels" begin
        @variables r θ
        @variables rs[1:nmodes] thetas[1:nmodes]
        rs_vec = collect(rs)
        thetas_vec = collect(thetas)
        n = rand(1:10)
        ns = rand(1:10, nmodes)
        test_configs = [(amplifier, (r, n), (rs_vec, ns), "Gaussian Amplifier"),(attenuator, (θ, n), (thetas_vec, ns), "Attenuator Channel")]
        for (op, single_params, multi_params, desc) in test_configs
            @testset "$desc" begin
                op_pair   = op(qpairbasis, single_params...)
                op_block  = op(qblockbasis, single_params...)
                @test op_pair isa GaussianChannel && op_block isa GaussianChannel
                @test op(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, single_params...) isa GaussianChannel
                @test op(Array, qpairbasis, single_params...) isa GaussianChannel
                @test isapprox(op(qblockbasis, single_params...), changebasis(QuadBlockBasis, op_pair))
            end
            @testset "$desc" begin
                op_pair_multi  = op(qpairbasis, multi_params...)
                op_block_multi = op(qblockbasis, multi_params...)
                @test op_pair_multi isa GaussianChannel && op_block_multi isa GaussianChannel
                @test op(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, multi_params...) isa GaussianChannel
                @test op(Array, qpairbasis, multi_params...) isa GaussianChannel
                @test isapprox(op(qblockbasis, multi_params...), changebasis(QuadBlockBasis, op_pair_multi))
            end
        end
    end
end

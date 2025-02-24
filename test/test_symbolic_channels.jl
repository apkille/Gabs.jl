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
        @variables alphas[1:2*nmodes]
        alphas_vec = collect(alphas)
        op_pair = displace(qpairbasis, α, noise)
        op_block = displace(qblockbasis, α, noise)
        @test op_pair isa GaussianChannel && op_block isa GaussianChannel
        @test displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, α, noise) isa GaussianChannel
        @test displace(Array, qpairbasis, α, noise) isa GaussianChannel
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
            end
            @testset "$desc" begin
                op_pair_multi  = op(qpairbasis, multi_params..., noise)
                op_block_multi = op(qblockbasis, multi_params..., noise)
                @test op_pair_multi isa GaussianChannel && op_block_multi isa GaussianChannel
                @test iszero(simplify(op_block_multi.transform - changebasis(QuadBlockBasis, op_pair_multi).transform))
                @test iszero(simplify(op_block_multi.disp - changebasis(QuadBlockBasis, op_pair_multi).disp))
                @test op(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, multi_params..., noise) isa GaussianChannel
                @test op(Array, qpairbasis, multi_params..., noise) isa GaussianChannel
            end
        end
    end

    @testset "Symbolic Gaussian Two-Mode Squeeze and Beamsplitter Operators" begin
        @variables r θ
        @variables rs[1:nmodes] thetas[1:nmodes]
        rs_vec = collect(rs)
        thetas_vec = collect(thetas)
        @variables θ_bs
        @variables thetas_bs[1:nmodes]
        thetas_bs_vec = collect(thetas_bs)
        test_configs = [(op = twosqueeze, single_params = (r, θ), multi_params = (rs_vec, thetas_vec), desc = "Two-Mode Squeeze Operator"),
                        (op = beamsplitter, single_params = (θ_bs,), multi_params = (thetas_bs_vec,), desc = "Beamsplitter Operator")]
        for (op, single_params, multi_params, desc) in test_configs
            @testset "$desc" begin
                op_single = op(2*qpairbasis, single_params..., noise_ds)
                @test op_single isa GaussianChannel
                @test op(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, single_params..., noise_ds) isa GaussianChannel
                @test op(Array, 2*qpairbasis, single_params..., noise_ds) isa GaussianChannel
                op_multi = op(2*qpairbasis, multi_params..., noise_ds)
                @test op_multi isa GaussianChannel
                @test op(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, multi_params..., noise_ds) isa GaussianChannel
                @test op(Array, 2*qpairbasis, multi_params..., noise_ds) isa GaussianChannel
            end
        end
    end

    @testset "Symbolic Gaussian Attenuator and Amplifier" begin
        @variables θ r
        @variables thetas[1:nmodes] rs[1:nmodes]
        n = rand(1:10)
        ns = rand(1:10, nmodes)
        for (op, param, param_vec) in [(attenuator, θ, thetas), (amplifier, r, rs)]
            op_pair = op(qpairbasis, param, n)
            op_block = op(qblockbasis, param, n)
            @test op_pair isa GaussianChannel && op_block isa GaussianChannel
            @test iszero(simplify(op_block.transform - changebasis(QuadBlockBasis, op_pair).transform))
            @test iszero(simplify(op_block.disp - changebasis(QuadBlockBasis, op_pair).disp))
            param_vec_col = collect(param_vec)
            ns_vec = collect(ns)
            op_pair_multi = op(qpairbasis, param_vec_col, ns_vec)
            op_block_multi = op(qblockbasis, param_vec_col, ns_vec)
            @test op_pair_multi isa GaussianChannel && op_block_multi isa GaussianChannel
            @test iszero(simplify(op_block_multi.transform - changebasis(QuadBlockBasis, op_pair_multi).transform))
            @test iszero(simplify(op_block_multi.disp - changebasis(QuadBlockBasis, op_pair_multi).disp))
        end
    end
end

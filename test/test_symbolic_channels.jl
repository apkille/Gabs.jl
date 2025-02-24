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

    @testset "Symbolic two-mode squeeze operator" begin
        @variables r θ
        @variables rs[1:nmodes] thetas[1:nmodes]
        rs_vec = collect(rs)
        thetas_vec = collect(thetas)
        op = twosqueeze(2*qpairbasis, r, θ, noise_ds)
        @test op isa GaussianChannel
        @test twosqueeze(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, r, θ, noise_ds) isa GaussianChannel
        @test twosqueeze(Array, 2*qpairbasis, r, θ, noise_ds) isa GaussianChannel
        @test isapprox(twosqueeze(2*qblockbasis, r, θ, noise_ds), changebasis(QuadBlockBasis, op))
        @test isapprox(twosqueeze(2*qblockbasis, rs_vec, thetas_vec, noise_ds), changebasis(QuadBlockBasis, twosqueeze(2*qpairbasis, rs_vec, thetas_vec, noise_ds)))
    end

    @testset "Symbolic beamsplitter operator" begin
        @variables θ_bs
        @variables thetas_bs[1:nmodes]
        thetas_bs_vec = collect(thetas_bs)
        op = beamsplitter(2*qpairbasis, theta, noise_ds)
        @test op isa GaussianChannel
        @test beamsplitter(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, θ_bs, noise_ds) isa GaussianChannel
        @test beamsplitter(Array, 2*qpairbasis, θ_bs, noise_ds) isa GaussianChannel
        @test isapprox(beamsplitter(2*qblockbasis, θ_bs, noise_ds), changebasis(QuadBlockBasis, op))
        @test isapprox(beamsplitter(2*qblockbasis, thetas_bs_vec, noise_ds), changebasis(QuadBlockBasis, beamsplitter(2*qpairbasis, thetas_bs_vec, noise_ds)))
    end

    @testset "Symbolic Gaussian Amplifier" begin
        @variables r
        @variables rs[1:nmodes]
        rs_vec = collect(rs)
        n = rand(1:10)
        ns = rand(1:10, nmodes)
        op_pair = amplifier(qpairbasis, r, n)
        op_block = amplifier(qblockbasis, r, n)
        @test op_pair isa GaussianChannel && op_block isa GaussianChannel
        @test amplifier(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, r, n) isa GaussianChannel
        @test amplifier(Array, qpairbasis, r, n) isa GaussianChannel
        @test isapprox(amplifier(qblockbasis, r, n), changebasis(QuadBlockBasis, op_pair))
        @test isapprox(amplifier(qblockbasis, rs_vec, ns), changebasis(QuadBlockBasis, amplifier(qpairbasis, rs_vec, ns)))
    end

    @testset "Symbolic attenuator channel" begin
        @variables θ
        @variables thetas[1:nmodes]
        thetas_vec = collect(thetas)
        n = rand(1:10)
        ns = rand(1:10, nmodes)
        op_pair = attenuator(qpairbasis, θ, n)
        op_block = attenuator(qblockbasis, θ, n)
        @test op_pair isa GaussianChannel && op_block isa GaussianChannel
        @test attenuator(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, θ, n) isa GaussianChannel
        @test attenuator(Array, qpairbasis, θ, n) isa GaussianChannel
        @test isapprox(attenuator(qblockbasis, θ, n), changebasis(QuadBlockBasis, op_pair))
        @test isapprox(attenuator(qblockbasis, thetas_vec, ns), changebasis(QuadBlockBasis, attenuator(qpairbasis, thetas_vec, ns)))
    end
end

@testitem "Symbolic Channels" begin
    using Gabs
    using Symbolics
    using StaticArrays

    nmodes = rand(1:5)
    qpairbasis = QuadPairBasis(nmodes)
    qblockbasis = QuadBlockBasis(nmodes)
    noise = rand(2*nmodes, 2*nmodes)
    noise_ds = [noise zeros(2*nmodes,2*nmodes); zeros(2*nmodes,2*nmodes) noise]
    T = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(2*nmodes), j in Base.OneTo(2*nmodes)
        if (j == 2*i-1) || (j + 2*nmodes == 2*i)
            T[i,j] = 1
        end
    end
    T_ds = zeros(4*nmodes, 4*nmodes)
    @inbounds for i in Base.OneTo(4*nmodes), j in Base.OneTo(4*nmodes)
        if (j == 2*i-1) || (j + 2*(2*nmodes) == 2*i)
            T_ds[i,j] = 1
        end
    end

    @testset "displacement operator" begin
        noise = zeros(Num, 2*nmodes, 2*nmodes)
        @variables α
        @variables alphas[1:2*nmodes]
        alphas_vec = collect(alphas)
        op_pair = displace(qpairbasis, α, noise)
        op_block = displace(qblockbasis, α, noise)
        @test op_pair isa GaussianChannel && op_block isa GaussianChannel
        @test displace(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, α, noise) isa GaussianChannel
        @test displace(Array, qpairbasis, α, noise) isa GaussianChannel
        @test isequal(displace(qblockbasis, α, T*noise*transpose(T)).disp, changebasis(QuadBlockBasis, op_pair).disp)
        @test isequal(displace(qblockbasis, α, T*noise*transpose(T)).basis, changebasis(QuadBlockBasis, op_pair).basis)
        @test isequal(displace(qblockbasis, α, T*noise*transpose(T)).transform, changebasis(QuadBlockBasis, op_pair).transform)
        @test isequal(displace(qblockbasis, alphas_vec, T*noise*transpose(T)).disp, changebasis(QuadBlockBasis, displace(qpairbasis, alphas_vec, noise)).disp)
        @test isequal(displace(qblockbasis, alphas_vec, T*noise*transpose(T)).basis, changebasis(QuadBlockBasis, displace(qpairbasis, alphas_vec, noise)).basis)
        @test isequal(displace(qblockbasis, alphas_vec, T*noise*transpose(T)).transform, changebasis(QuadBlockBasis, displace(qpairbasis, alphas_vec, noise)).transform)
    end

    @testset "Symbolic Squeeze and Phase-Shift Operators" begin
        noise = zeros(Num, 2*nmodes, 2*nmodes)
        @variables r θs
        @variables rs[1:nmodes] thetas_s[1:nmodes]
        rs_vec = collect(rs)
        thetas_vec = collect(thetas_s)
        @variables θps
        @variables thetas_ps[1:nmodes]
        thetas_ps_vec = collect(thetas_ps)
        sub_dict(params_tuple) = Dict(vcat([p .=> 0 for p in params_tuple]...)...)
        test_configs = [(squeeze, (r, θs), (rs_vec, thetas_vec), "Squeeze Operator"), (phaseshift, (θps,), (thetas_ps_vec,), "Phase-Shift Operator")]
        for (op, single_params, multi_params, desc) in test_configs
            @testset "Single mode: $desc" begin
                op_pair   = op(qpairbasis, single_params..., noise)
                op_block  = op(qblockbasis, single_params..., noise)
                @test op_pair isa GaussianChannel && op_block isa GaussianChannel
                @test iszero(simplify(substitute(op_block.transform - changebasis(QuadBlockBasis, op_pair).transform, sub_dict(single_params))))
                @test iszero(simplify(op_block.disp - changebasis(QuadBlockBasis, op_pair).disp))
                @test op(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, single_params..., noise) isa GaussianChannel
                @test op(Array, qpairbasis, single_params..., noise) isa GaussianChannel
                @test isequal(op(qblockbasis, single_params..., T*noise*transpose(T)).disp, changebasis(QuadBlockBasis, op_pair).disp)
                @test isequal(op(qblockbasis, single_params..., T*noise*transpose(T)).basis, changebasis(QuadBlockBasis, op_pair).basis)
                @test isequal(op(qblockbasis, single_params..., T*noise*transpose(T)).transform, changebasis(QuadBlockBasis, op_pair).transform)
          end
            @testset "Multi mode: $desc" begin
                op_pair_multi  = op(qpairbasis, multi_params..., noise)
                op_block_multi = op(qblockbasis, multi_params..., noise)
                @test op_pair_multi isa GaussianChannel && op_block_multi isa GaussianChannel
                @test iszero(simplify(substitute(op_block_multi.transform - changebasis(QuadBlockBasis, op_pair_multi).transform, sub_dict(multi_params))))
                @test iszero(simplify(op_block_multi.disp - changebasis(QuadBlockBasis, op_pair_multi).disp))
                @test op(SVector{2*nmodes}, SMatrix{2*nmodes,2*nmodes}, qpairbasis, multi_params..., noise) isa GaussianChannel
                @test op(Array, qpairbasis, multi_params..., noise) isa GaussianChannel
                @test isequal(op(qblockbasis, multi_params..., T*noise*transpose(T)).disp, changebasis(QuadBlockBasis, op(qpairbasis, multi_params..., noise)).disp)
                @test isequal(op(qblockbasis, multi_params..., T*noise*transpose(T)).basis, changebasis(QuadBlockBasis, op(qpairbasis, multi_params..., noise)).basis)
                @test isequal(op(qblockbasis, multi_params..., T*noise*transpose(T)).transform, changebasis(QuadBlockBasis, op(qpairbasis, multi_params..., noise)).transform)
            end
        end
    end

    @testset "Symbolic Gaussian Two-Mode Squeeze  and Beamsplitter Operators" begin
        noise = zeros(Num, 2*nmodes, 2*nmodes)
        noise_ds = [noise zeros(Num, 2*nmodes, 2*nmodes); zeros(Num, 2*nmodes, 2*nmodes) noise]
        @variables r θ
        @variables rs[1:nmodes] thetas[1:nmodes]
        rs_vec = collect(rs)
        thetas_vec = collect(thetas)
        @variables θ_bs
        @variables thetas_bs[1:nmodes]
        thetas_bs_vec = collect(thetas_bs)
        sub_dict(params_tuple) = Dict(vcat([p .=> 0 for p in params_tuple]...)...)
        test_configs = [(op = twosqueeze, single_params = (r, θ), multi_params = (rs_vec, thetas_vec), desc = "Two-Mode Squeeze Operator"),
                        (op = beamsplitter, single_params = (θ_bs,), multi_params = (thetas_bs_vec,), desc = "Beamsplitter Operator")]
        for (op, single_params, multi_params, desc) in test_configs
            @testset "$desc" begin
                # Test single-mode operation
                op_single = op(2*qpairbasis, single_params..., noise_ds)
                @test op_single isa GaussianChannel
                @test iszero(simplify(substitute(op_single.transform - changebasis(QuadBlockBasis, op_single).transform, sub_dict(single_params))))
                @test iszero(simplify(op_single.disp - changebasis(QuadBlockBasis, op_single).disp))
                @test op(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, single_params..., noise_ds) isa GaussianChannel
                @test op(Array, 2*qpairbasis, single_params..., noise_ds) isa GaussianChannel
                @test isequal(op(2*qblockbasis, single_params..., T_ds*noise_ds*transpose(T_ds)).disp, changebasis(QuadBlockBasis, op_single).disp)
                @test isequal(op(2*qblockbasis, single_params..., T_ds*noise_ds*transpose(T_ds)).basis, changebasis(QuadBlockBasis, op_single).basis)
                @test isequal(op(2*qblockbasis, single_params..., T_ds*noise_ds*transpose(T_ds)).transform, changebasis(QuadBlockBasis, op_single).transform)
                # Test multi-mode operation
                op_multi = op(2*qpairbasis, multi_params..., noise_ds)
                @test op_multi isa GaussianChannel
                @test iszero(simplify(substitute(op_multi.transform - changebasis(QuadBlockBasis, op_multi).transform, sub_dict(multi_params))))
                @test iszero(simplify(op_multi.disp - changebasis(QuadBlockBasis, op_multi).disp))
                @test op(SVector{4*nmodes}, SMatrix{4*nmodes,4*nmodes}, 2*qpairbasis, multi_params..., noise_ds) isa GaussianChannel
                @test op(Array, 2*qpairbasis, multi_params..., noise_ds) isa GaussianChannel
                @test isequal(op(2*qblockbasis, multi_params..., T_ds*noise_ds*transpose(T_ds)).disp, changebasis(QuadBlockBasis, op_multi).disp)
                @test isequal(op(2*qblockbasis, multi_params..., T_ds*noise_ds*transpose(T_ds)).basis, changebasis(QuadBlockBasis, op_multi).basis)
                @test isequal(op(2*qblockbasis, multi_params..., T_ds*noise_ds*transpose(T_ds)).transform, changebasis(QuadBlockBasis, op_multi).transform)
            end
        end
    end

    @testset "Symbolic Gaussian Attenuator and Amplifier" begin
        @variables θ r
        @variables thetas[1:nmodes] rs[1:nmodes]
        n = rand(1:10)
        ns = rand(1:10, nmodes)
        noise_ds = [noise zeros(Num, 2*nmodes,2*nmodes); zeros(Num, 2*nmodes,2*nmodes) noise]
        for (op, param, param_vec) in [(attenuator, θ, thetas), (amplifier, r, rs)]
            # Test single mode
            op_pair = op(qpairbasis, param, n)
            op_block = op(qblockbasis, param, n)
            @test op_pair isa GaussianChannel && op_block isa GaussianChannel
            @test iszero(simplify(substitute(op_block.transform - changebasis(QuadBlockBasis, op_pair).transform, Dict(vcat(param.=> 0)...))))
            @test iszero(simplify(op_block.disp - changebasis(QuadBlockBasis, op_pair).disp))
            # Test multi-mode
            param_vec_col = collect(param_vec)
            ns_vec = collect(ns)
            op_pair_multi = op(qpairbasis, param_vec_col, ns_vec)
            op_block_multi = op(qblockbasis, param_vec_col, ns_vec)
            @test op_pair_multi isa GaussianChannel && op_block_multi isa GaussianChannel
            @test iszero(simplify(substitute(op_block_multi.transform - changebasis(QuadBlockBasis, op_pair_multi).transform, Dict(vcat(param_vec_col.=> 0)...))))
            @test iszero(simplify(op_block_multi.disp - changebasis(QuadBlockBasis, op_pair_multi).disp))
        end
    end
end

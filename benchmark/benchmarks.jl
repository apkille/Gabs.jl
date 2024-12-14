using BenchmarkTools
using Gabs
using LinearAlgebra
using PkgBenchmark


const SUITE = BenchmarkGroup()

prob_list = ("operations")
op_list = ("unitary product", "channel product", "tensor product", "partial trace")
for prob in prob_list
    SUITE[prob] = BenchmarkGroup([prob])
end
for op in op_list
    SUITE["operations"][op] = BenchmarkGroup()
end

for nmodes in (2, 10, 50, 100, 200)
    basis = QuadPairBasis(nmodes)
    SUITE["operations"]["unitary product"][string(nmodes)] = @benchmarkable op * state setup=(op = randunitary($basis); state = randstate($basis))
    SUITE["operations"]["channel product"][string(nmodes)] = @benchmarkable ch * state setup=(ch = randchannel($basis); state = randstate($basis))
    SUITE["operations"]["tensor product"][string(nmodes)] = @benchmarkable state1 ⊗ state2 setup=(state1 = randstate($basis); state2 = randstate($basis))
    SUITE["operations"]["partial trace"][string(nmodes)] = @benchmarkable ptrace(state, indices) setup=(state = randstate($basis); indices = collect(2:$nmodes))
end
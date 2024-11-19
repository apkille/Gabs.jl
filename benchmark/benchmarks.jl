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
    SUITE["operations"]["unitary product"][string(nmodes)] = @benchmarkable op * state setup=(op = randunitary($nmodes); state = randstate($nmodes))
    SUITE["operations"]["channel product"][string(nmodes)] = @benchmarkable ch * state setup=(ch = randchannel($nmodes); state = randstate($nmodes))
    SUITE["operations"]["tensor product"][string(nmodes)] = @benchmarkable state1 ⊗ state2 setup=(state1 = randstate($nmodes); state2 = randstate($nmodes))
    SUITE["operations"]["partial trace"][string(nmodes)] = @benchmarkable ptrace(state, indices) setup=(state = randstate($nmodes); indices = collect(2:$nmodes))
end
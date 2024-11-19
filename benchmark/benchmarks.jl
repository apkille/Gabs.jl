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

for dim in (2, 10, 50, 100, 200)
    SUITE["operations"]["unitary product"][string(dim)] = @benchmarkable op * state setup=(op = randunitary($dim); state = randstate($dim))
    SUITE["operations"]["channel product"][string(dim)] = @benchmarkable ch * state setup=(ch = randchannel($dim); state = randstate($dim))
    SUITE["operations"]["tensor product"][string(dim)] = @benchmarkable state1 ⊗ state2 setup=(state1 = randstate($dim); state2 = randstate($dim))
    SUITE["operations"]["partial trace"][string(dim)] = @benchmarkable ptrace(state, collect(2:dim)) setup=(state = randstate($dim))
end
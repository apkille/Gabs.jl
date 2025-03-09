module StaticArraysExt

using StaticArrays: SVector, SMatrix, SArray

using Gabs
using Gabs: SymplecticBasis, QuadPairBasis,_randchannel

import Gabs: ptrace, tensor, ⊗, _promote_output_matrix, _promote_output_vector,
            _generaldyne_map, _infer_types, randchannel

include("utils.jl")
include("measurements.jl")

end

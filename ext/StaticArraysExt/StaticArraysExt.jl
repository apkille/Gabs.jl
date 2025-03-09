module StaticArraysExt

using StaticArrays: SVector, SMatrix, SArray

using Gabs
using Gabs: SymplecticBasis, QuadPairBasis

import Gabs: ptrace, tensor, ⊗, _promote_output_matrix, _promote_output_vector,
            _generaldyne_map, _infer_types

include("utils.jl")
include("measurements.jl")

end

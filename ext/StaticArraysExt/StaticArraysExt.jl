module StaticArraysExt

using StaticArrays: SVector, SMatrix, SArray

using Gabs
using Gabs: SymplecticBasis, QuadPairBasis

import Gabs: ptrace, tensor, ⊗, _promote_output_matrix, _promote_output_vector,
            _generaldyne_map, infer_mean_type, infer_covar_type, infer_displacement_type,
            infer_symplectic_type

include("utils.jl")
include("measurements.jl")

end

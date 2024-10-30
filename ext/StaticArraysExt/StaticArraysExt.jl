module StaticArraysExt

using StaticArrays: SVector, SMatrix, SArray

using Gabs
import Gabs: ptrace, tensor, ⊗, _promote_output_matrix, _promote_output_vector, 
            _generaldyne_map

include("utils.jl")
include("measurements.jl")

end
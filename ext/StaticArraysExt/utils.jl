Base.@propagate_inbounds function _promote_output_vector(::Type{T1}, ::Type{T2}, vec_out) where {T1<:SVector, T2<:SVector}
    lengths = length(T1) + length(T2)
    return SVector{lengths}(vec_out)
end
Base.@propagate_inbounds function _promote_output_vector(::Type{<:SVector}, vec_out, vec_length::Int)
    return SVector{vec_length}(vec_out)
end
Base.@propagate_inbounds function _promote_output_matrix(::Type{T1}, ::Type{T2}, mat_out) where {T1<:SMatrix,T2<:SMatrix}
    dim = size(T1)[1] + size(T2)[1]
    return SMatrix{dim,dim}(mat_out)
end
Base.@propagate_inbounds function _promote_output_matrix(::Type{<:SMatrix}, mat_out, out_dim::Int)
    return SMatrix{out_dim,out_dim}(mat_out)
end
Base.@propagate_inbounds function _promote_output_matrix(::Type{<:SMatrix}, mat_out, out_dim::Tuple)
    return SMatrix{out_dim[1],out_dim[2]}(mat_out)
end

function infer_mean_type(::Type{SVector}, basis::Gabs.SymplecticBasis{N}) where {N}
    nmodes = basis.nmodes
    return SVector{2*nmodes, Float64}
end

function infer_covar_type(::Type{SMatrix}, basis::Gabs.SymplecticBasis{N}) where {N}
    nmodes = basis.nmodes
    return SMatrix{2*nmodes, 2*nmodes, Float64, 4*nmodes*nmodes}
end

function infer_displacement_type(::Type{SVector}, basis::Gabs.QuadPairBasis{N}) where {N}
    nmodes = basis.nmodes
    return SVector{2*nmodes, Float64}
end

function infer_symplectic_type(::Type{SMatrix}, basis::Gabs.QuadPairBasis{N}) where {N}
    nmodes = basis.nmodes
    return SMatrix{2*nmodes, 2*nmodes, Float64, 4*nmodes*nmodes}
end

function infer_transform_type(::Type{SMatrix}, basis::Gabs.SymplecticBasis{N}) where {N}
    nmodes = basis.nmodes
    return SMatrix{2*nmodes, 2*nmodes, Float64, 4*nmodes*nmodes}
end

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
Base.@propagate_inbounds function _promote_output_matrix(::Type{T1}, ::Type{T2}, mat_out) where {T1<:SMatrix, T2<:AbstractMatrix}
    return SMatrix{size(mat_out,1), size(mat_out,2)}(mat_out)
end
Base.@propagate_inbounds function _promote_output_vector(::Type{T1}, ::Type{T2}, vec_out) where {T1<:SVector, T2<:AbstractVector}
    return SVector{length(vec_out)}(vec_out)
end

function _infer_types(::Type{T1}, ::Type{T2}, basis) where {T1<:SVector, T2<:SMatrix}
    nmodes = basis.nmodes
    eltypeT1, eltypeT2 = eltype(T1), eltype(T2)
    return SVector{2*nmodes, eltypeT1}, SMatrix{2*nmodes, 2*nmodes, eltypeT2}
end
function _infer_types(::Type{T1}, ::Type{T2}, basis) where {T1<:SVector{<:Int}, T2<:SMatrix{<:Int, <:Int}}
    nmodes = basis.nmodes
    return T1, T2
end
function _infer_types(::Type{T}, basis) where {T<:SArray}
    nmodes = basis.nmodes
    eltypeT = eltype(T)
    return SVector{2*nmodes, eltypeT}, SMatrix{2*nmodes, 2*nmodes, eltypeT}
end

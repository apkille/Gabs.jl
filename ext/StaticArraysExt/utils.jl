function _promote_output_vector(::Type{T1}, ::Type{T2}, vec_out) where {T1<:SArray, T2<:SArray}
    lengths = length(T1) + length(T2)
    return SVector{lengths}(vec_out)
end
function _promote_output_vector(::Type{T}, vec_out, vec_length::Tl) where {T<:SArray,Tl<:Int}
    return SVector{vec_length}(vec_out)
end
function _promote_output_matrix(::Type{T1}, ::Type{T2}, mat_out) where {T1<:SArray,T2<:SArray}
    dim = size(T1)[1] + size(T2)[1]
    return SMatrix{dim,dim}(mat_out)
end
function _promote_output_matrix(::Type{T}, mat_out, out_dim::Td) where {T<:SArray,Td<:Int}
    return SMatrix{out_dim, out_dim}(mat_out)
end
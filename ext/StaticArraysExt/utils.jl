function _promote_output_vector(::Type{T1}, ::Type{T2}, vec_out) where {T1<:SArray,T2}
    T = promote_type(T1, T2)
    if T <: SArray
        lengths = length(T1) + length(T2)
        vec_out′ = SVector{lengths}(vec_out)
    elseif T <: Vector
        vec_out′ = vec_out
    else
        vec_out′ = T(vec_out)
    end
    return vec_out′
end
function _promote_output_vector(::Type{T1}, ::Type{T2}, vec_out) where {T1,T2<:SArray}
    T = promote_type(T1, T2)
    if T <: SArray
        lengths = length(T1) + length(T2)
        vec_out′ = SVector{lengths}(vec_out)
    elseif T <: Vector
        vec_out′ = vec_out
    else
        vec_out′ = T(vec_out)
    end
    return vec_out′
end
function _promote_output_vector(::Type{T1}, ::Type{T2}, vec_out) where {T1<:SArray, T2<:SArray}
    lengths = length(T1) + length(T2)
    return SVector{lengths}(vec_out)
end
function _promote_output_vector(::Type{T}, vec_out, vec_length::Tl) where {T<:SArray,Tl<:Int}
    return SVector{vec_length}(vec_out)
end
function _promote_output_matrix(::Type{T1}, ::Type{T2}, mat_out) where {T1<:SArray,T2}
    T = promote_type(T1, T2)
    if T <: SArray
        dim = size(T1)[1] + size(T2)[1]
        mat_out′ = SMatrix{dim,dim}(mat_out)
    elseif T <: Matrix
        mat_out′ = mat_out
    else
        mat_out′ = T(mat_out)
    end
    return mat_out′
end
function _promote_output_matrix(::Type{T1}, ::Type{T2}, mat_out) where {T1,T2<:SArray}
    T = promote_type(T1, T2)
    if T <: SArray
        dim = size(T1)[1] + size(T2)[1]
        mat_out′ = SMatrix{dim,dim}(mat_out)
    elseif T <: Matrix
        mat_out′ = mat_out
    else
        mat_out′ = T(mat_out)
    end
    return mat_out′
end
function _promote_output_matrix(::Type{T1}, ::Type{T2}, mat_out) where {T1<:SArray,T2<:SArray}
    dim = size(T1)[1] + size(T2)[1]
    return SMatrix{dim,dim}(mat_out)
end
function _promote_output_matrix(::Type{T}, mat_out, out_dim::Td) where {T<:SArray,Td<:Int}
    return SMatrix{out_dim, out_dim}(mat_out)
end
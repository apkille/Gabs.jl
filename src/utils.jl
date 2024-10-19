function _promote_output_vector(::Type{T1}, ::Type{T2}, vec_out) where {T1,T2}
    T = promote_type(T1, T2)
    if T <: Vector
        vec_out′ = vec_out
    else
        vec_out′ = T(vec_out)
    end
    return vec_out′
end
function _promote_output_vector(::Type{T}, vec_out, vec_length::Tl) where {T,Tl<:Int}
    if T <: Vector
        vec_out′ = vec_out
    else
        vec_out′ = T(vec_out)
    end
    return vec_out′
end
function _promote_output_matrix(::Type{T1}, ::Type{T2}, mat_out) where {T1,T2}
    T = promote_type(T1, T2)
    if T <: Matrix
        mat_out′ = mat_out
    else
        mat_out′ = T(mat_out)
    end
    return mat_out′
end
function _promote_output_matrix(::Type{T}, mat_out, out_dim::Td) where {T,Td<:Int}
    if T <: Matrix
        mat_out′ = mat_out
    else
        mat_out′ = T(mat_out)
    end
    return mat_out′
end
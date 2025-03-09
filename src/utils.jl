Base.@propagate_inbounds function _promote_output_vector(::Type{T1}, ::Type{T2}, vec_out) where {T1,T2}
    T = promote_type(T1, T2)
    T <: Vector{Float64} ? vec_out : T(vec_out)
end
Base.@propagate_inbounds function _promote_output_vector(::Type{T}, vec_out, vec_length::Tl) where {T,Tl<:Int}
    T <: Vector{Float64} ? vec_out : T(vec_out)
end
Base.@propagate_inbounds function _promote_output_matrix(::Type{T1}, ::Type{T2}, mat_out) where {T1,T2}
    T = promote_type(T1, T2)
    T <: Matrix{Float64} ? mat_out : T(mat_out)
end
Base.@propagate_inbounds function _promote_output_matrix(::Type{T}, mat_out, out_dim::Td) where {T,Td<:Int}
    T <: Matrix{Float64} ? mat_out : T(mat_out)
end
Base.@propagate_inbounds function _promote_output_matrix(::Type{T}, mat_out, out_dim::Td) where {T,Td<:Tuple}
    T <: Matrix{Float64} ? mat_out : T(mat_out)
end

function _infer_types end

_infer_types(T1, T2, basis) = T1
_infer_types(T, basis) = T

function _infer_types(::Type{Array{T}}, basis) where {T}
    nmodes = basis.nmodes
    return (Vector{T}, Matrix{T})
end
function _infer_types(::Type{Array}, basis)
    return _infer_types(Array{Float64}, basis)
end
function _infer_types(::Type{Vector{T}}, ::Type{Matrix{T}}, basis) where {T}
    nmodes = basis.nmodes
    return (Vector{T}, Matrix{T})
end
function _infer_types(::Type{Vector}, ::Type{Matrix}, basis)
    return _infer_types(Vector{Float64}, Matrix{Float64}, basis)
end

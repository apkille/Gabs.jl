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
    return collect(vec_out)
end

Base.@propagate_inbounds function _promote_output_vector(::Type{T1}, ::Type{T2}, vec_out) where {T1<:AbstractVector, T2<:SVector}
    return collect(vec_out)
end

Base.@propagate_inbounds function _promote_output_vector(::Type{Vector}, vec_out, vec_length::Int)
    return Vector{eltype(vec_out)}(vec_out)
end

abstract type ArrayTrait end

struct DenseArrayTrait <: ArrayTrait end
struct StaticArrayTrait <: ArrayTrait end

array_trait(::Type{<:Array}) = DenseArrayTrait()
array_trait(::Type{<:SArray}) = StaticArrayTrait()
array_trait(::Type{<:UnionAll}) = StaticArrayTrait()

function _infer_types(::DenseArrayTrait, nmodes, T = Float64)
    disp_type = Vector{T}
    transform_type = Matrix{T}
    return disp_type, transform_type
end
function _infer_types(::StaticArrayTrait, nmodes, T = Float64)
    disp_type = SArray{Tuple{2*nmodes}, T}
    transform_type = SArray{Tuple{2*nmodes, 2*nmodes}, T}
    return disp_type, transform_type
end
function _infer_types(T1, T2, basis)
    nmodes = basis.nmodes
    elT1 = eltype(T1)
    elT2 = eltype(T2)
    disp_type1, _ = _infer_types(array_trait(T1), nmodes, elT1)
    _, transform_type2 = _infer_types(array_trait(T2), nmodes, elT2)
    return disp_type1, transform_type2
end
function _infer_types(T, basis)
    nmodes = basis.nmodes
    elT = eltype(T)
    disp_type, transform_type = _infer_types(array_trait(T), nmodes, elT)
    return disp_type, transform_type
end

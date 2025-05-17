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

# % random.jl

randunitary(::Type{SArray}, basis::SymplecticBasis{N}; passive = false, ħ = 2) where {N<:Int} = 
    let n = basis.nmodes
        randunitary(SVector{2n,Float64}, SMatrix{2n,2n,Float64}, basis; passive = passive, ħ = ħ)
    end
randunitary(::Type{SVector}, ::Type{SMatrix}, basis::SymplecticBasis{N}; passive = false, ħ = 2) where {N<:Int} = 
    let n = basis.nmodes
        randunitary(SVector{2n,Float64}, SMatrix{2n,2n,Float64}, basis; passive = passive, ħ = ħ)
    end
randunitary(::Type{SVector{M}}, ::Type{SMatrix{M,M}}, basis::SymplecticBasis{N}; passive = false, ħ = 2) where {M,N<:Int} = 
    let n = basis.nmodes
        M == 2n || error("Size mismatch: SVector{$M}/SMatrix{$M,$M} != 2n (n=$n)")
        randunitary(SVector{2n,Float64}, SMatrix{2n,2n,Float64}, basis; passive = passive, ħ = ħ)
    end

randchannel(::Type{SArray}, basis::SymplecticBasis{N}) where {N<:Int} =
    let n = basis.nmodes
        randchannel(SVector{2n,Float64}, SMatrix{2n,2n,Float64}, basis)
    end
randchannel(::Type{SVector}, ::Type{SMatrix}, basis::SymplecticBasis{N}) where {N<:Int} =
    let n = basis.nmodes
        randchannel(SVector{2n,Float64}, SMatrix{2n,2n,Float64}, basis)
    end
randchannel(::Type{SVector{M}}, ::Type{SMatrix{M,M}}, basis::SymplecticBasis{N}) where {M,N<:Int} =
    let n = basis.nmodes
        M == 2n || error("Size mismatch: SVector{$M}/SMatrix{$M,$M} != 2n (n=$n)")
        randchannel(SVector{2n,Float64}, SMatrix{2n,2n,Float64}, basis)
    end

randsymplectic(::Type{SArray}, basis::SymplecticBasis{N}; passive = false) where {N<:Int} = 
    let n = basis.nmodes
        randsymplectic(SMatrix{2n,2n,Float64}, basis; passive = passive)
    end
randsymplectic(::Type{SMatrix}, basis::SymplecticBasis{N}; passive = false) where {N<:Int} =
    let n = basis.nmodes
        randsymplectic(SMatrix{2n,2n,Float64}, basis; passive = passive)
    end
randsymplectic(::Type{SMatrix{M,M}}, basis::SymplecticBasis{N}; passive = false) where {M,N<:Int} =
    let n = basis.nmodes
        M == 2n || error("Size mismatch: SMatrix{$M,$M} != 2n×2n (n=$n)")
        randsymplectic(SMatrix{2n,2n,Float64}, basis; passive = passive)
    end

randstate(::Type{SArray}, basis::SymplecticBasis{N}; pure = false, ħ = 2) where {N<:Int} =
    let n = basis.nmodes
        randstate(SVector{2n,Float64}, SMatrix{2n,2n,Float64}, basis; pure = pure, ħ = ħ)
    end
randstate(::Type{SVector}, ::Type{SMatrix}, basis::SymplecticBasis{N}; pure = false, ħ = 2) where {N<:Int} =
    let n = basis.nmodes
        randstate(SVector{2n,Float64}, SMatrix{2n,2n,Float64}, basis; pure = pure, ħ = ħ)
    end
randstate(::Type{SVector{M}}, ::Type{SMatrix{M,M}}, basis::SymplecticBasis{N}; pure = false, ħ = 2) where {M,N<:Int} =
    let n = basis.nmodes
        M == 2n || error("Size mismatch: SVector{$M}/SMatrix{$M,$M} != 2n (n=$n)")
        randstate(SVector{2n,Float64}, SMatrix{2n,2n,Float64}, basis; pure = pure, ħ = ħ)
    end

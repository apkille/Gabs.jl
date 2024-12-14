abstract type SymplecticBasis{N} end

"""
Defines a symplectic basis for a bosonic system of size `nmodes` in which
the quadrature field operators are arranged pairwise.
"""
struct QuadPairBasis{N} <: SymplecticBasis{N}
    nmodes::N
end

"""
Defines a symplectic basis for a bosonic system of size `nmodes` in which
the quadrature field operators are arranged blockwise.
"""
struct QuadBlockBasis{N} <: SymplecticBasis{N}
    nmodes::N
end

function Base.:(*)(n::N, repr2::R) where {N<:Number,R<:SymplecticBasis}
    R(n*repr2.nmodes)
end
function Base.:(+)(repr1::R, repr2::R) where {R<:SymplecticBasis}
    R(repr1.nmodes + repr2.nmodes)
end
function Base.:(-)(repr1::R, repr2::R) where {R<:SymplecticBasis}
    R(repr1.nmodes - repr2.nmodes)
end

"""
    symplecticform([T = Matrix{Float64},] basis::SymplecticBasis)

Compute the symplectic form matrix of size 2N x 2N corresponding to `basis`.
"""
function symplecticform(basis::QuadPairBasis{N}) where {N<:Int}
    nmodes = basis.nmodes
    Omega = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        Omega[2*i-1, 2*i] = 1.0
        Omega[2*i, 2*i-1] = -1.0
    end
    return Omega
end
function symplecticform(basis::QuadBlockBasis{N}) where {N<:Int}
    nmodes = basis.nmodes
    Omega = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in 1:nmodes, j in nmodes:2*nmodes
        if isequal(i, j-nmodes)
            Omega[i,j] = 1.0
        end
    end
    @inbounds for i in nmodes:2*nmodes, j in 1:nmodes
        if isequal(i-nmodes,j)
            Omega[i, j] = -1.0
        end
    end
    return Omega
end
symplecticform(::Type{T}, basis::SymplecticBasis{N}) where {T, N<:Int} = T(symplecticform(basis))
"""
    williamson(state::GaussianState) -> Williamson

Compute the williamson decomposition of the `covar` field
of a `Gaussian state` object and return a `Williamson` object.

A symplectic matrix `S` and symplectic spectrum `spectrum` can be obtained
via `F.S` and `F.spectrum`.

Iterating the decomposition produces the components `S` and `spectrum`.

To compute only the symplectic spectrum of a Gaussian state, call [`sympspectrum`](@ref).
"""
function williamson(x::GaussianState{<:QuadBlockBasis,M,V}) where {M,V}
    basis = x.basis
    return williamson(BlockForm(basis.nmodes), x.covar)
end
function williamson(x::GaussianState{<:QuadPairBasis,M,V}) where {M,V}
    basis = x.basis
    return williamson(PairForm(basis.nmodes), x.covar)
end

"""
    polar(state::GaussianState) -> Polar

Compute the Polar decomposition of the `symplectic` field
of a `Gaussian unitary` object and return a `Polar` object.

`O` and `P` can be obtained from the factorization `F` via `F.O` and `F.P`, such that `S = O * P`.
For the symplectic polar decomposition case, `O` is an orthogonal symplectic matrix and `P` is a positive-definite
symmetric symplectic matrix.

Iterating the decomposition produces the components `O` and `P`.
"""
function polar(x::GaussianUnitary{<:QuadBlockBasis,D,S}) where {D,S}
    basis = x.basis
    return polar(x.symplectic)
end
function polar(x::GaussianUnitary{<:QuadPairBasis,D,S}) where {D,S}
    basis = x.basis
    return polar(x.symplectic)
end
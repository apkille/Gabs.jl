"""
    williamson(state::GaussianState) -> Williamson

Compute the williamson decomposition of the `covar` field
of a `GaussianState` object and return a `Williamson` object.

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
    polar(state::GaussianUnitary) -> Polar

Compute the Polar decomposition of the `symplectic` field
of a `GaussianUnitary` object and return a `Polar` object.

`O` and `P` can be obtained from the factorization `F` via `F.O` and `F.P`, such that `S = O * P`.
For the symplectic polar decomposition case, `O` is an orthogonal symplectic matrix and `P` is a positive-definite
symmetric symplectic matrix.

Iterating the decomposition produces the components `O` and `P`.
"""
polar(x::GaussianUnitary{B,D,S}) where {B<:SymplecticBasis,D,S} = polar(x.symplectic)

"""
    blochmessiah(state::GaussianUnitary) -> BlochMessiah

Compute the Bloch-Messiah/Euler decomposition of the `symplectic` field
of a `GaussianUnitary` and return a `BlockMessiah` object.

The orthogonal symplectic matrices `O` and `Q` as well as the singular values `values` can be obtained
via `F.O`, `F.Q`, and `F.values`, respectively.

Iterating the decomposition produces the components `O`, `values`, and `Q`, in that order.
"""
function blochmessiah(x::GaussianUnitary{<:QuadBlockBasis,D,S}) where {D,S}
    basis = x.basis
    return blochmessiah(BlockForm(basis.nmodes), x.symplectic)
end
function blochmessiah(x::GaussianUnitary{<:QuadPairBasis,D,S}) where {D,S}
    basis = x.basis
    return blochmessiah(PairForm(basis.nmodes), x.symplectic)
end
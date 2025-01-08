"""
    williamson(state::GaussianState) -> Williamson

Compute the williamson decomposition of the `covar` field
of a `Gaussian state` object and return a `Williamson` object.

A symplectic matrix `S` and symplectic spectrum `spectrum` can be obtained
via `F.S` and `F.spectrum`.

Iterating the decomposition produces the components `S` and `spectrum`.
"""
function williamson(x::GaussianState{<:QuadBlockBasis,M,V}) where {M,V}
    basis = x.basis
    return williamson(BlockForm(basis.nmodes), x.covar)
end
function williamson(x::GaussianState{<:QuadPairBasis,M,V}) where {M,V}
    basis = x.basis
    return williamson(PairForm(basis.nmodes), x.covar)
end